import torch as th 
import triton
import triton.language as tl

# wrapper for cuda code, allow using optimized CUDA kernel easily from python. 
from torch.utils.cpp_extension import load
import os 
package_dir = os.path.dirname(os.path.abspath(__file__))
sseqs_sw_ext = load(name="sseqs_sw_ext", sources=[f"{package_dir}/sw_bind.cpp", f"{package_dir}/sw.cu"], extra_cuda_cflags=["-O3", "--use_fast_math"],  extra_cflags=["-O3"])
def sw(query: str, targets: list[str], gap_open=11, gap_extend=1):
    _q_tensor = th.tensor([AA_MAP[aa] for aa in query if aa in AA_MAP], dtype=th.uint8, device="cuda")
    ascii = th.hstack([th.tensor([ord(c) for c in target+'@'], dtype=th.uint8, device="cuda") for target in targets])
    delimiter = 64
    starts = th.nonzero(ascii == delimiter, as_tuple=False).flatten()
    good_idx = th.arange(len(targets), dtype=th.int32, device="cuda")
    return sseqs_sw_ext.sw_cuda_affine(
        _q_tensor,
        good_idx,
        ascii,
        starts.to(th.int32),
        gap_open=gap_open,
        gap_extend=gap_extend
    )

# @alex: code below used to build .a3m file not search DB.
# it needs backtracking -- this is just quick&dirty triton, should rewrite to CUDA. 
def sw_affine_backtrack(query: str, targets: list[str], gap_open=11, gap_extend=1, device="cuda", scores=False):
    """
    Computes Smith-Waterman alignment scores with Affine Gap Penalties using BLOSUM62.
    Returns M_matrix_full, best_scores_tensor, and aligned_sequences_list.
    Gap penalties should be positive values representing costs (e.g., gap_open=11, gap_extend=1).
    """
    B = len(targets)
    Qrow = encode_seq(query, device).unsqueeze(0)
    Q_tensor = Qrow.expand(B, -1).contiguous()
    T_tensor, t_lens, W_max_target = pack_targets(targets, device)

    query_actual_len = Q_tensor.shape[1] - 1
    Q_COLS_dim = Q_tensor.shape[1]
    T_COLS_dim = T_tensor.shape[1]

    buffer_dim = query_actual_len + 1
    WIDTH = 1 << (buffer_dim - 1).bit_length() if buffer_dim > 0 else 1
    if query_actual_len == 0: WIDTH = max(1, WIDTH)

    prev2_M = th.zeros((B, WIDTH), dtype=th.int16, device=device)
    prev1_M = th.zeros_like(prev2_M)
    curr_M  = th.zeros_like(prev2_M)

    prev2_Ix = th.zeros_like(prev2_M)
    prev1_Ix = th.zeros_like(prev2_M)
    curr_Ix  = th.zeros_like(prev2_M)

    prev2_Iy = th.zeros_like(prev2_M)
    prev1_Iy = th.zeros_like(prev2_M)
    curr_Iy  = th.zeros_like(prev2_M)
    
    best_scores_tensor = th.zeros(B, dtype=th.int16, device=device) # Stores max H scores
    th.cuda.empty_cache()
    
    # so i need to keep all of this in memory? 
    M_matrix_full  = th.zeros((B, Q_COLS_dim, T_COLS_dim), dtype=th.int16, device=device)
    Ix_matrix_full = th.zeros_like(M_matrix_full)
    Iy_matrix_full = th.zeros_like(M_matrix_full)
    #ic('sw_affine_backtrack', M_matrix_full.nbytes/1e9, Ix_matrix_full.nbytes/1e9, Iy_matrix_full.nbytes/1e9)
    
    t_actual_lens_vec = t_lens - 1

    max_actual_target_len_for_dmax = 0
    if t_actual_lens_vec.numel() > 0:
        valid_target_lengths_for_dmax = t_actual_lens_vec[t_actual_lens_vec >= 0]
        if valid_target_lengths_for_dmax.numel() > 0:
            max_actual_target_len_for_dmax = int(valid_target_lengths_for_dmax.max().item())
    
    if query_actual_len == 0 and max_actual_target_len_for_dmax == 0:
        d_max_val = 2
    else:
        d_max_val = query_actual_len + max_actual_target_len_for_dmax + 1

    blosum_matrix_tensor = th.tensor(_BLOSUM62_FLAT_LIST, dtype=th.int32, device=device)
    char_to_index_map_tensor = th.tensor(_CHAR_TO_BLOSUM_IDX_LIST, dtype=th.int32, device=device)
    AA_SIZE_CONST = len(_AA_ORDER)
    # DEFAULT_SUB_PENALTY_CONST is defined globally

    # before for loop

    for d_val in range(2, d_max_val): # can we have this be
        sw_diag_batch_affine[(B,)]( 
            prev2_M, prev1_M,
            prev1_Ix, prev2_Ix, 
            prev1_Iy, prev2_Iy,
            Q_tensor, T_tensor,
            d_val,
            gap_open, gap_extend, 
            blosum_matrix_tensor, char_to_index_map_tensor,
            curr_M, curr_Ix, curr_Iy, 
            best_scores_tensor, t_lens,
            M_matrix_full, Ix_matrix_full, Iy_matrix_full, 
            Q_COLS=Q_COLS_dim, # Use keyword argument
            T_COLS=T_COLS_dim, # Use keyword argument
            Q_LEN_ACTUAL=query_actual_len,
            WIDTH=WIDTH,
            AA_SIZE=AA_SIZE_CONST,
            DEFAULT_SUB_PENALTY=DEFAULT_SUB_PENALTY_CONST
        )

        temp_p1_M = prev1_M
        prev1_M = curr_M
        prev2_M = temp_p1_M
        curr_M = th.zeros_like(prev1_M, device=device) 

        temp_p1_Ix = prev1_Ix
        prev1_Ix = curr_Ix
        prev2_Ix = temp_p1_Ix
        curr_Ix = th.zeros_like(prev1_Ix, device=device)

        temp_p1_Iy = prev1_Iy
        prev1_Iy = curr_Iy
        prev2_Iy = temp_p1_Iy
        curr_Iy = th.zeros_like(prev1_Iy, device=device)

    # forloop

    max_target_len_in_batch = 0
    if t_actual_lens_vec.numel() > 0:
        valid_lengths = t_actual_lens_vec[t_actual_lens_vec >= 0]
        if valid_lengths.numel() > 0:
            max_target_len_in_batch = valid_lengths.max().item()
    
    MAX_ALIGN_LEN = query_actual_len + max_target_len_in_batch + 1 
    if MAX_ALIGN_LEN == 0: MAX_ALIGN_LEN = 1 

    out_q_aligned_chars = th.zeros((B, MAX_ALIGN_LEN), dtype=th.int32, device=device)
    out_t_aligned_chars = th.zeros((B, MAX_ALIGN_LEN), dtype=th.int16, device=device)
    out_q_aligned_len   = th.zeros(B, dtype=th.int32, device=device)
    out_t_aligned_len   = th.zeros(B, dtype=th.int32, device=device)
    out_q_start_0based  = th.zeros(B, dtype=th.int32, device=device)
    out_q_end_0based    = th.zeros(B, dtype=th.int32, device=device)
    out_t_start_0based  = th.zeros(B, dtype=th.int32, device=device)
    out_t_end_0based    = th.zeros(B, dtype=th.int32, device=device)

    GAP_CHAR_CODE_CONST = ord('-')
    _MAX_ORD_VAL_PYTHON = 128 
    VERY_NEGATIVE_SCORE_CONST_PY = -16384 

    M_matrix_full = M_matrix_full.to(th.int32)  # Ensure M_matrix_full is int32 for consistency
    Ix_matrix_full = Ix_matrix_full.to(th.int32)
    Iy_matrix_full = Iy_matrix_full.to(th.int32)

    backtrack_kernel_batched_affine[(B,)](
        M_matrix_full, Ix_matrix_full, Iy_matrix_full,
        Q_tensor, T_tensor, best_scores_tensor, t_actual_lens_vec,
        query_actual_len, 
        gap_open, gap_extend, 
        Q_COLS_dim, T_COLS_dim,
        blosum_matrix_tensor, char_to_index_map_tensor,
        AA_SIZE_CONST, DEFAULT_SUB_PENALTY_CONST,
        out_q_aligned_chars, out_t_aligned_chars,
        out_q_aligned_len, out_t_aligned_len,
        out_q_start_0based, out_q_end_0based,
        out_t_start_0based, out_t_end_0based,
        MAX_ALIGN_LEN=MAX_ALIGN_LEN,
        GAP_CHAR_CODE=GAP_CHAR_CODE_CONST,
        _MAX_ORD_VAL_CONSTEXPR=_MAX_ORD_VAL_PYTHON,
        VERY_NEGATIVE_SCORE_CONST=VERY_NEGATIVE_SCORE_CONST_PY
    ) # presumably we can half this with int16? or just make one kernel for this? 
    #backtrack kernel batched affine

    out_q_aligned_chars_cpu = out_q_aligned_chars.cpu()
    out_t_aligned_chars_cpu = out_t_aligned_chars.cpu()
    out_q_aligned_len_cpu = out_q_aligned_len.cpu()

    # 1.  Bulk-clean every byte on the C side (no Python branch per element)
    #     Anything <0 or >255 becomes ord('?') == 63
    q_buf = th.where(
        (out_q_aligned_chars_cpu < 0) | (out_q_aligned_chars_cpu > 255),
        th.full_like(out_q_aligned_chars_cpu, 63),      # '?'
        out_q_aligned_chars_cpu,
    ).to(th.uint8).contiguous().numpy()                 # (B, MAX_ALIGN_LEN)  uint8

    t_buf = th.where(
        (out_t_aligned_chars_cpu < 0) | (out_t_aligned_chars_cpu > 255),
        th.full_like(out_t_aligned_chars_cpu, 63),
        out_t_aligned_chars_cpu,
    ).to(th.uint8).contiguous().numpy()

    # 2.  Bring all the scalar columns across once
    q_len_arr   = out_q_aligned_len_cpu.cpu().numpy()
    q_start_arr = out_q_start_0based.cpu().numpy()
    q_end_arr   = out_q_end_0based.cpu().numpy()
    t_start_arr = out_t_start_0based.cpu().numpy()
    t_end_arr   = out_t_end_0based.cpu().numpy()
    scores_arr  = best_scores_tensor.cpu().numpy()

    # 3.  Fast construction ---------------------------------------------------------
    aligned_sequences_list_fast = []
    for b in range(B):
        q_len = q_len_arr[b]
        if q_len:
            s = MAX_ALIGN_LEN - q_len
            # NumPy view → bytes view → str; zero copy until the very last step
            q_aligned_str = q_buf[b, s : s + q_len].tobytes().decode("latin-1")
            t_aligned_str = t_buf[b, s : s + q_len].tobytes().decode("latin-1")
        else:
            q_aligned_str = t_aligned_str = ""

        aligned_sequences_list_fast.append(
            {
                "q_aligned": q_aligned_str,
                "t_aligned": t_aligned_str,
                "q_start_orig_0based": int(q_start_arr[b]),
                "q_end_orig_0based":   int(q_end_arr[b]),
                "t_start_orig_0based": int(t_start_arr[b]),
                "t_end_orig_0based":   int(t_end_arr[b]),
                "score": float(scores_arr[b]),
            }
        )
    if scores: return best_scores_tensor.to(th.int32)
    return M_matrix_full, best_scores_tensor.to(th.int32), aligned_sequences_list_fast




# Define BLOSUM62 matrix and AA mapping (can be global or passed to sw10)
AA_MAP = {aa: idx for idx, aa in enumerate("ARNDCQEGHILKMFPSTWYV")}
_AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
_AA_MAP_ORD_IDX = {aa: idx for idx, aa in enumerate(_AA_ORDER)}

# BLOSUM62 matrix values 
_BLOSUM62_FLAT_LIST = [
    #A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
     4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0, # A
    -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3, # R
    -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, # N
    -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, # D
     0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1, # C
    -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, # Q
    -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, # E
     0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3, # G
    -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, # H
    -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3, # I
    -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1, # L
    -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, # K
    -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1, # M
    -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1, # F
    -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2, # P
     1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, # S
     0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0, # T
    -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3, # W
    -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1, # Y
     0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4  # V
]

_MAX_ORD_VAL = 128 # Sufficient for ASCII characters
_CHAR_TO_BLOSUM_IDX_LIST = [-1] * _MAX_ORD_VAL
for char, idx in _AA_MAP_ORD_IDX.items():
    if ord(char) < _MAX_ORD_VAL:
        _CHAR_TO_BLOSUM_IDX_LIST[ord(char)] = idx

DEFAULT_SUB_PENALTY_CONST = -5 # Ensure this is defined globally or passed appropriately

def encode_seq(s: str, device="cuda") -> th.IntTensor:
    v = th.zeros(len(s) + 1, dtype=th.int16, device=device)   # leading 0
    v[1:] = th.tensor(list(map(ord, s)), dtype=th.int16)
    return v

def pack_targets(targets: list[str], device="cuda"):
    W = max(map(len, targets)) + 1
    B = len(targets)
    T = th.zeros((B, W), dtype=th.int16, device=device)
    lens = th.empty(B, dtype=th.int32,  device=device)
    for i, t in enumerate(targets):
        T[i, 1:len(t)+1] = th.tensor(list(map(ord, t)), dtype=th.int16)
        lens[i] = len(t) + 1
    return T, lens, W

@triton.jit
def sw_diag_batch(prev2, prev1, Q, T,
                  d,
                  gap,
                  blosum_matrix_ptr, char_to_index_map_ptr,
                  out, best_scores, tlen_ptr,
                  DP_matrix_ptr,
                  Q_COLS: tl.constexpr, T_COLS: tl.constexpr,
                  Q_LEN_ACTUAL: tl.constexpr,
                  WIDTH: tl.constexpr,
                  AA_SIZE: tl.constexpr,
                  DEFAULT_SUB_PENALTY: tl.constexpr):

    pid   = tl.program_id(0)                     # which target in batch
    tlen  = tl.load(tlen_ptr  + pid) # len(target_seq_k) + 1
    t_actual_len = tlen - 1 # M_k

    # --- Calculate i_min and m internally ---
    # i_min_k = max(1, d - M_k)
    # i_max_k = min(N, d - 1)
    # m_k = i_max_k - i_min_k + 1
    # Ensure d-1 is not negative for i_max calculation
    d_minus_1 = d - 1
    if d_minus_1 < 0: # Should not happen with d starting at 2
        d_minus_1 = 0

    i_min_val = d - t_actual_len
    if i_min_val < 1:
        i_min_val = 1
    i_min = i_min_val

    i_max_val = Q_LEN_ACTUAL # N
    if i_max_val > d_minus_1:
        i_max_val = d_minus_1
    i_max = i_max_val
    
    m = i_max - i_min + 1
    # --- End of i_min and m calculation ---

    offs  = tl.arange(0, WIDTH)
    # active_cell_mask: only process if m > 0 and offs < m
    active_cell_mask = (m > 0) & (offs < m)

    # Pointers for diagonal buffers (prev2, prev1, out)
    P2_base = prev2 + pid * WIDTH
    P1_base = prev1 + pid * WIDTH
    Out_base = out  + pid * WIDTH

    # Pointers to the start of sequence data for Q and T for this batch item
    # Q is (B, Q_COLS), T is (B, T_COLS)
    # Q_COLS = len(query) + 1
    # T_COLS = max_target_len + 1 (overall padded width for T)
    Q_batch_base = Q + pid * Q_COLS
    T_batch_base = T + pid * T_COLS

    # Load scores from previous diagonals
    diag_score_from_h = tl.load(P2_base + (i_min - 1) + offs, active_cell_mask, other=0)
    up_score_from_h   = tl.load(P1_base + (i_min - 1) + offs, active_cell_mask, other=0)
    left_score_from_h = tl.load(P1_base + (i_min    ) + offs, active_cell_mask, other=0)

    # Current DP cell coordinates (1-based for matrix, and for 1-padded sequences)
    i_dp = i_min + offs # Query sequence dimension index (1 to len(query))
    j_dp = d - i_dp    # Target sequence dimension index (1 to len(target_k))

    # --- Character loading ---
    # For Q: Q[0] is pad, Q[1] is 1st char. Q_COLS = len(query) + 1. Valid char indices in Q: 1 to Q_COLS-1.
    q_char_load_mask = active_cell_mask & (i_dp >= 1) & (i_dp < Q_COLS)
    ai = tl.load(Q_batch_base + i_dp, q_char_load_mask, other=0) # Load Q[i_dp] (ord(char) or 0 for pad)

    # For T: T[batch_idx, 0] is pad. tlen = len(target_k) + 1. Valid char indices in T for this target: 1 to tlen-1.
    t_char_load_mask = active_cell_mask & (j_dp >= 1) & (j_dp < tlen)
    bj = tl.load(T_batch_base + j_dp, t_char_load_mask, other=0) # Load T[j_dp] (ord(char) or 0 for pad)

    # --- Substitution score using BLOSUM62 ---
    # Map ord(char) to BLOSUM indices (0 to AA_SIZE-1, or -1 if not an AA/padding)
    # char_to_index_map_ptr should map ord(0) (padding) to -1.
    idx_ai = tl.load(char_to_index_map_ptr + ai, active_cell_mask, other=-1)
    idx_bj = tl.load(char_to_index_map_ptr + bj, active_cell_mask, other=-1)

    # Determine if both characters are valid AAs for BLOSUM lookup
    blosum_lookup_mask = (idx_ai != -1) & (idx_bj != -1) & active_cell_mask
    
    # Load score from BLOSUM matrix if valid AA pair, otherwise use default penalty
    sub_score_val = tl.load(blosum_matrix_ptr + idx_ai * AA_SIZE + idx_bj,
                            mask=blosum_lookup_mask,
                            other=DEFAULT_SUB_PENALTY).to(tl.int16)
    
    current_score = tl.maximum(tl.zeros_like(diag_score_from_h), 
                               tl.maximum(diag_score_from_h + sub_score_val,
                                          tl.maximum(up_score_from_h + gap, left_score_from_h + gap)))

    # Store the computed score for the current cell into the 'out' buffer (for next diagonal calculation)
    tl.store(Out_base + i_min + offs, current_score, active_cell_mask)

    # --- Store to full DP matrix M[pid, i_dp, j_dp] ---
    # i_dp and j_dp are 1-based for DP matrix cells.
    # Q_LEN_ACTUAL is N (actual query length). t_actual_len is M_k (actual target length for this item).
    # DP matrix M is (B, Q_COLS, T_COLS) where Q_COLS = N+1, T_COLS = Max_M_k+1.
    # Valid indices for M[pid, :, :] are M[pid, 1..N, 1..M_k].
    
    dp_store_mask = active_cell_mask & \
                    (i_dp >= 1) & (i_dp <= Q_LEN_ACTUAL) & \
                    (j_dp >= 1) & (j_dp <= t_actual_len) # Use t_actual_len for this specific target

    if m > 0: # Only proceed if there are cells on this diagonal for this batch item
        # Calculate flat offset into DP_matrix_ptr
        # DP_matrix_ptr is dimensioned (B, Q_COLS, T_COLS)
        dp_offset = pid * Q_COLS * T_COLS + \
                    i_dp * T_COLS + \
                    j_dp
        tl.store(DP_matrix_ptr + dp_offset, current_score, mask=dp_store_mask)
    # --- End of DP matrix store ---
    
    # Find the maximum score in the current diagonal for this batch item
    # Use a block reduction to find the maximum
    current_max = tl.max(tl.where(active_cell_mask, current_score, tl.zeros_like(current_score)))
    
    # Update the best score for this batch item
    best_score = tl.load(best_scores + pid)
    best_score = tl.maximum(best_score, current_max)
    tl.store(best_scores + pid, best_score)
    
    # Copy out to prev2 for next iteration (prev1 becomes prev2, out becomes prev1)
    # We don't need to zero out 'out' buffer here as we'll overwrite values in the next iteration
    for i in range(0, WIDTH):
        if i < WIDTH:  # Always true, but helps triton optimize
            tl.store(P2_base + i, tl.load(P1_base + i))
            tl.store(P1_base + i, tl.load(Out_base + i))
            tl.store(Out_base + i, 0)  # Zero out for next iteration

# ────────────────── AFFINE GAP SW KERNEL ───────────────────────
@triton.jit
def sw_diag_batch_affine(
    # Diagonal buffers for M (match/mismatch) scores
    prev2_M_ptr, prev1_M_ptr,
    # Diagonal buffers for Ix (insertion in Q) scores
    prev1_Ix_ptr, prev2_Ix_ptr, # prev2_Ix for M_ij = S + max(M_i-1,j-1, Ix_i-1,j-1, Iy_i-1,j-1)
    # Diagonal buffers for Iy (insertion in T) scores
    prev1_Iy_ptr, prev2_Iy_ptr, # prev2_Iy for M_ij = S + max(M_i-1,j-1, Ix_i-1,j-1, Iy_i-1,j-1)
    # Sequences
    Q_ptr, T_ptr,
    # Current diagonal index
    d,
    # Gap penalties
    gap_open, gap_extend,
    # Scoring data
    blosum_matrix_ptr, char_to_index_map_ptr,
    # Output diagonal buffers
    out_M_ptr, out_Ix_ptr, out_Iy_ptr,
    # Overall best scores per item
    best_scores_ptr,
    # Target sequence lengths (including padding char)
    tlen_ptr,
    # Full DP matrices for M, Ix, and Iy scores
    M_full_matrix_ptr, Ix_full_matrix_ptr, Iy_full_matrix_ptr,
    # Dimensions and constants
    Q_COLS: tl.constexpr, T_COLS: tl.constexpr, # Padded lengths: N+1, Max_M_k+1
    Q_LEN_ACTUAL: tl.constexpr, # Actual query length N
    WIDTH: tl.constexpr,        # Width of diagonal buffers
    AA_SIZE: tl.constexpr,
    DEFAULT_SUB_PENALTY: tl.constexpr
):
    pid = tl.program_id(0)
    tlen = tl.load(tlen_ptr + pid)
    t_actual_len = tlen - 1

    d_minus_1 = d - 1
    if d_minus_1 < 0: d_minus_1 = 0
    i_min_val = d - t_actual_len
    if i_min_val < 1: i_min_val = 1
    i_min = i_min_val
    i_max_val = Q_LEN_ACTUAL
    if i_max_val > d_minus_1: i_max_val = d_minus_1
    i_max = i_max_val
    m = i_max - i_min + 1

    offs = tl.arange(0, WIDTH)
    active_cell_mask = (m > 0) & (offs < m)

    # Base pointers for diagonal buffers for this batch item
    P2_M_base  = prev2_M_ptr  + pid * WIDTH
    P1_M_base  = prev1_M_ptr  + pid * WIDTH
    Out_M_base = out_M_ptr    + pid * WIDTH

    P2_Ix_base  = prev2_Ix_ptr + pid * WIDTH
    P1_Ix_base  = prev1_Ix_ptr + pid * WIDTH
    Out_Ix_base = out_Ix_ptr   + pid * WIDTH

    P2_Iy_base  = prev2_Iy_ptr + pid * WIDTH
    P1_Iy_base  = prev1_Iy_ptr + pid * WIDTH
    Out_Iy_base = out_Iy_ptr   + pid * WIDTH
    
    Q_batch_base = Q_ptr + pid * Q_COLS
    T_batch_base = T_ptr + pid * T_COLS

    # Load predecessor scores for M(i,j) calculation
    # These are M(i-1,j-1), Ix(i-1,j-1), Iy(i-1,j-1) from diagonal d-2
    m_val_d2  = tl.load(P2_M_base  + (i_min - 1) + offs, active_cell_mask, other=0)
    ix_val_d2 = tl.load(P2_Ix_base + (i_min - 1) + offs, active_cell_mask, other=0)
    iy_val_d2 = tl.load(P2_Iy_base + (i_min - 1) + offs, active_cell_mask, other=0)

    # Load predecessor scores for Ix(i,j) calculation
    # These are M(i-1,j), Ix(i-1,j), Iy(i-1,j) from diagonal d-1
    m_val_d1_for_ix  = tl.load(P1_M_base  + (i_min - 1) + offs, active_cell_mask, other=0) # M[i-1,j]
    ix_val_d1_for_ix = tl.load(P1_Ix_base + (i_min - 1) + offs, active_cell_mask, other=0) # Ix[i-1,j]
    iy_val_d1_for_ix = tl.load(P1_Iy_base + (i_min - 1) + offs, active_cell_mask, other=0) # Iy[i-1,j]

    # Load predecessor scores for Iy(i,j) calculation
    # These are M(i,j-1), Ix(i,j-1), Iy(i,j-1) from diagonal d-1
    m_val_d1_for_iy  = tl.load(P1_M_base  + i_min + offs, active_cell_mask, other=0) # M[i,j-1]
    ix_val_d1_for_iy = tl.load(P1_Ix_base + i_min + offs, active_cell_mask, other=0) # Ix[i,j-1]
    iy_val_d1_for_iy = tl.load(P1_Iy_base + i_min + offs, active_cell_mask, other=0) # Iy[i,j-1]

    i_dp = i_min + offs
    j_dp = d - i_dp

    q_char_load_mask = active_cell_mask & (i_dp >= 1) & (i_dp < Q_COLS)
    ai = tl.load(Q_batch_base + i_dp, q_char_load_mask, other=0)
    t_char_load_mask = active_cell_mask & (j_dp >= 1) & (j_dp < tlen)
    bj = tl.load(T_batch_base + j_dp, t_char_load_mask, other=0)

    idx_ai = tl.load(char_to_index_map_ptr + ai, active_cell_mask, other=-1)
    idx_bj = tl.load(char_to_index_map_ptr + bj, active_cell_mask, other=-1)
    blosum_lookup_mask = (idx_ai != -1) & (idx_bj != -1) & active_cell_mask
    sub_score_val = tl.load(blosum_matrix_ptr + idx_ai * AA_SIZE + idx_bj,
                            mask=blosum_lookup_mask,
                            other=DEFAULT_SUB_PENALTY).to(tl.int16)

    # Calculate M(i,j) = s(i,j) + H(i-1,j-1)
    # H(i-1,j-1) = max(0, M(i-1,j-1), Ix(i-1,j-1), Iy(i-1,j-1))
    h_val_prev_diag = tl.maximum(tl.zeros_like(m_val_d2),
                                 tl.maximum(m_val_d2, tl.maximum(ix_val_d2, iy_val_d2)))
    current_M_val = sub_score_val + h_val_prev_diag
    
    # Calculate Ix(i,j) = max( H(i-1,j) - gap_open, Ix(i-1,j) - gap_extend )
    # H(i-1,j) = max(0, M(i-1,j), Ix(i-1,j), Iy(i-1,j))
    h_val_up = tl.maximum(tl.zeros_like(m_val_d1_for_ix),
                          tl.maximum(m_val_d1_for_ix, tl.maximum(ix_val_d1_for_ix, iy_val_d1_for_ix)))
    score_ix_from_h = h_val_up - gap_open

    VERY_NEGATIVE_SCORE = tl.full((), -16384, dtype=tl.int16)
    # If extending from Ix(0,j), treat Ix(0,j) as -infinity
    # i_dp is 1-based current query index. So i_dp-1 is query index of cell (i-1,j)
    is_ix_source_on_q_border = ((i_dp - 1) == 0)
    effective_ix_pred_for_extend = tl.where(is_ix_source_on_q_border, VERY_NEGATIVE_SCORE, ix_val_d1_for_ix)
    score_ix_from_ix = effective_ix_pred_for_extend - gap_extend
    current_Ix_val = tl.maximum(score_ix_from_h, score_ix_from_ix)

    # Calculate Iy(i,j) = max( H(i,j-1) - gap_open, Iy(i,j-1) - gap_extend )
    # H(i,j-1) = max(0, M(i,j-1), Ix(i,j-1), Iy(i,j-1))
    h_val_left = tl.maximum(tl.zeros_like(m_val_d1_for_iy),
                           tl.maximum(m_val_d1_for_iy, tl.maximum(ix_val_d1_for_iy, iy_val_d1_for_iy)))
    score_iy_from_h = h_val_left - gap_open

    # If extending from Iy(i,0), treat Iy(i,0) as -infinity
    # j_dp is 1-based current target index. So j_dp-1 is target index of cell (i,j-1)
    is_iy_source_on_t_border = ((j_dp - 1) == 0)
    effective_iy_pred_for_extend = tl.where(is_iy_source_on_t_border, VERY_NEGATIVE_SCORE, iy_val_d1_for_iy)
    score_iy_from_iy = effective_iy_pred_for_extend - gap_extend
    current_Iy_val = tl.maximum(score_iy_from_h, score_iy_from_iy)
    
    # Store M, Ix, Iy to output buffers for this diagonal
    tl.store(Out_M_base + i_min + offs, current_M_val, active_cell_mask)
    tl.store(Out_Ix_base + i_min + offs, current_Ix_val, active_cell_mask)
    tl.store(Out_Iy_base + i_min + offs, current_Iy_val, active_cell_mask)

    # Calculate final Smith-Waterman score H(i,j) = max(0, M(i,j), Ix(i,j), Iy(i,j))
    current_H_score = tl.maximum(tl.zeros_like(current_M_val),
                                 tl.maximum(current_M_val,
                                            tl.maximum(current_Ix_val, current_Iy_val)))
    
    # Store M(i,j), Ix(i,j), Iy(i,j) to full DP matrices
    dp_store_mask = active_cell_mask & \
                    (i_dp >= 1) & (i_dp <= Q_LEN_ACTUAL) & \
                    (j_dp >= 1) & (j_dp <= t_actual_len)
    if m > 0:
        dp_offset = pid * Q_COLS * T_COLS + i_dp * T_COLS + j_dp
        # The following three lines replace the previous single store to DP_matrix_ptr
        tl.store(M_full_matrix_ptr  + dp_offset, current_M_val, mask=dp_store_mask)
        tl.store(Ix_full_matrix_ptr + dp_offset, current_Ix_val, mask=dp_store_mask)
        tl.store(Iy_full_matrix_ptr + dp_offset, current_Iy_val, mask=dp_store_mask)

    # Update best score for this batch item (based on H score)
    current_max_H_on_diag = tl.max(tl.where(active_cell_mask, current_H_score, tl.zeros_like(current_H_score)))
    best_score_old = tl.load(best_scores_ptr + pid)
    best_score_new = tl.maximum(best_score_old, current_max_H_on_diag)
    tl.store(best_scores_ptr + pid, best_score_new)





@triton.jit
def backtrack_kernel_batched_affine(
    # Input DP matrices & sequences
    M_matrix_ptr, Ix_matrix_ptr, Iy_matrix_ptr, # (B, Q_COLS, T_COLS)
    Q_ptr, T_ptr,                               # (B, Q_COLS), (B, T_COLS)
    best_scores_ptr,                            # (B) - Best H score for each item
    t_actual_lens_ptr,                          # (B) - Actual length of each target sequence (M_k)
    query_actual_len,                           # scalar, N (actual query length)
    gap_open, gap_extend,                       # scalar
    Q_COLS_dim, T_COLS_dim,                     # scalar, N+1, Max_M_k+1

    # BLOSUM data
    blosum_matrix_ptr, char_to_index_map_ptr,
    AA_SIZE: tl.constexpr,
    DEFAULT_SUB_PENALTY: tl.constexpr,

    # Output buffers
    out_q_aligned_chars_ptr, out_t_aligned_chars_ptr,
    out_q_aligned_len_ptr, out_t_aligned_len_ptr,
    out_q_start_0based_ptr, out_q_end_0based_ptr,
    out_t_start_0based_ptr, out_t_end_0based_ptr,

    # Constants
    MAX_ALIGN_LEN: tl.constexpr,
    GAP_CHAR_CODE: tl.constexpr,
    _MAX_ORD_VAL_CONSTEXPR: tl.constexpr,
    VERY_NEGATIVE_SCORE_CONST: tl.constexpr
):
    pid = tl.program_id(0)

    # Initialize outputs for this batch item (same as linear kernel)
    tl.store(out_q_aligned_len_ptr + pid, 0)
    tl.store(out_t_aligned_len_ptr + pid, 0)
    tl.store(out_q_start_0based_ptr + pid, -1)
    tl.store(out_q_end_0based_ptr + pid, -1)
    tl.store(out_t_start_0based_ptr + pid, -1)
    tl.store(out_t_end_0based_ptr + pid, -1)
    
    item_best_score = tl.load(best_scores_ptr + pid)
    item_target_actual_len = tl.load(t_actual_lens_ptr + pid)
    if item_target_actual_len < 0:
        item_target_actual_len = 0

    can_align = ((item_best_score > 0) and (query_actual_len > 0)) and (item_target_actual_len > 0)

    if can_align:
        # --- 1. Find end cell (q_idx_end_1based, t_idx_end_1based) of local alignment ---
        q_idx_end_1based = -1
        t_idx_end_1based = -1
        found_end_cell_flag = False

        for i_loop_idx in range(query_actual_len): # 0 to N-1
            if not found_end_cell_flag:
                for j_loop_idx in range(item_target_actual_len): # 0 to M_k-1
                    if not found_end_cell_flag:
                        i_1based = i_loop_idx + 1
                        j_1based = j_loop_idx + 1

                        m_offset  = pid * Q_COLS_dim * T_COLS_dim + i_1based * T_COLS_dim + j_1based
                        val_m     = tl.load(M_matrix_ptr  + m_offset)
                        val_ix    = tl.load(Ix_matrix_ptr + m_offset)
                        val_iy    = tl.load(Iy_matrix_ptr + m_offset)
                        
                        val_h = tl.maximum(0, tl.maximum(val_m, tl.maximum(val_ix, val_iy)))

                        if val_h == item_best_score:
                            q_idx_end_1based = i_1based
                            t_idx_end_1based = j_1based
                            found_end_cell_flag = True
        
        if found_end_cell_flag: 
            Q_item_ptr = Q_ptr + pid * Q_COLS_dim
            T_item_ptr = T_ptr + pid * T_COLS_dim
            
            curr_q_char_write_idx = MAX_ALIGN_LEN - 1
            curr_t_char_write_idx = MAX_ALIGN_LEN - 1
            current_alignment_len = 0

            curr_i_1based = q_idx_end_1based
            curr_j_1based = t_idx_end_1based
            
            keep_tracing_flag = True
            for _trace_iter in range(query_actual_len + item_target_actual_len + 2) :
                if keep_tracing_flag: 
                    current_offset = pid * Q_COLS_dim * T_COLS_dim + curr_i_1based * T_COLS_dim + curr_j_1based
                    current_m_val  = tl.load(M_matrix_ptr + current_offset)
                    current_ix_val = tl.load(Ix_matrix_ptr + current_offset)
                    current_iy_val = tl.load(Iy_matrix_ptr + current_offset)
                    current_h_val  = tl.maximum(0, tl.maximum(current_m_val, tl.maximum(current_ix_val, current_iy_val)))

                    if ((current_h_val == 0) or (curr_i_1based <= 0)) or (curr_j_1based <= 0): 
                        keep_tracing_flag = False
                        # Removed continue, outer loop's 'if keep_tracing_flag' will prevent further processing this iteration
                    
                    # Nested check: only proceed with main traceback logic if flag is still true
                    if keep_tracing_flag:
                        ord_char_q = tl.load(Q_item_ptr + curr_i_1based)
                        ord_char_t = tl.load(T_item_ptr + curr_j_1based)
                        
                        moved_this_step = False

                        if current_h_val == current_m_val:
                            mask_ord_q_valid = (ord_char_q >= 0) & (ord_char_q < _MAX_ORD_VAL_CONSTEXPR)
                            idx_q_char = tl.load(char_to_index_map_ptr + ord_char_q, mask=mask_ord_q_valid, other=-1)
                            mask_ord_t_valid = (ord_char_t >= 0) & (ord_char_t < _MAX_ORD_VAL_CONSTEXPR)
                            idx_t_char = tl.load(char_to_index_map_ptr + ord_char_t, mask=mask_ord_t_valid, other=-1)
                            blosum_load_mask = (idx_q_char != -1) & (idx_t_char != -1)
                            s_ij = tl.load(blosum_matrix_ptr + idx_q_char * AA_SIZE + idx_t_char,
                                           mask=blosum_load_mask, other=DEFAULT_SUB_PENALTY).to(tl.int32)

                            h_prev_diag = tl.zeros_like(current_m_val) 
                            if (curr_i_1based > 1) and curr_j_1based > 1: 
                                prev_diag_offset = pid*Q_COLS_dim*T_COLS_dim + (curr_i_1based-1)*T_COLS_dim + (curr_j_1based-1)
                                m_prev_diag  = tl.load(M_matrix_ptr  + prev_diag_offset)
                                ix_prev_diag = tl.load(Ix_matrix_ptr + prev_diag_offset)
                                iy_prev_diag = tl.load(Iy_matrix_ptr + prev_diag_offset)
                                h_prev_diag  = tl.maximum(tl.full((), 0, dtype=tl.int32), tl.maximum(m_prev_diag, tl.maximum(ix_prev_diag, iy_prev_diag)))
                            
                            if current_m_val == s_ij + h_prev_diag:
                                tl.store(out_q_aligned_chars_ptr + pid * MAX_ALIGN_LEN + curr_q_char_write_idx, ord_char_q)
                                tl.store(out_t_aligned_chars_ptr + pid * MAX_ALIGN_LEN + curr_t_char_write_idx, ord_char_t)
                                curr_i_1based -= 1
                                curr_j_1based -= 1
                                moved_this_step = True

                        elif not moved_this_step and current_h_val == current_ix_val:
                            h_up = tl.zeros_like(current_ix_val) 
                            if curr_i_1based > 1 and curr_j_1based >= 1: 
                                up_offset = pid*Q_COLS_dim*T_COLS_dim + (curr_i_1based-1)*T_COLS_dim + curr_j_1based
                                m_up  = tl.load(M_matrix_ptr  + up_offset)
                                ix_up = tl.load(Ix_matrix_ptr + up_offset)
                                iy_up = tl.load(Iy_matrix_ptr + up_offset)
                                h_up  = tl.maximum(tl.full((), 0, dtype=tl.int32), tl.maximum(m_up, tl.maximum(ix_up, iy_up)))
                            
                            score_ix_from_h_up = h_up - gap_open
                            
                            ix_prev_up = tl.full((), VERY_NEGATIVE_SCORE_CONST, dtype=tl.int32)
                            if curr_i_1based > 1 and curr_j_1based >= 1:
                                ix_prev_up_offset = pid*Q_COLS_dim*T_COLS_dim + (curr_i_1based-1)*T_COLS_dim + curr_j_1based
                                ix_prev_up = tl.load(Ix_matrix_ptr + ix_prev_up_offset)

                            score_ix_from_ix_up = ix_prev_up - gap_extend

                            if current_ix_val == score_ix_from_h_up or \
                               (current_ix_val == score_ix_from_ix_up and score_ix_from_h_up < score_ix_from_ix_up) : 
                                tl.store(out_q_aligned_chars_ptr + pid * MAX_ALIGN_LEN + curr_q_char_write_idx, ord_char_q)
                                tl.store(out_t_aligned_chars_ptr + pid * MAX_ALIGN_LEN + curr_t_char_write_idx, GAP_CHAR_CODE)
                                curr_i_1based -= 1
                                moved_this_step = True
                            elif current_ix_val == score_ix_from_ix_up:
                                tl.store(out_q_aligned_chars_ptr + pid * MAX_ALIGN_LEN + curr_q_char_write_idx, ord_char_q)
                                tl.store(out_t_aligned_chars_ptr + pid * MAX_ALIGN_LEN + curr_t_char_write_idx, GAP_CHAR_CODE)
                                curr_i_1based -= 1
                                moved_this_step = True
                        
                        elif not moved_this_step and current_h_val == current_iy_val:
                            h_left = tl.zeros_like(current_iy_val)
                            if curr_j_1based > 1 and curr_i_1based >= 1:
                                left_offset = pid*Q_COLS_dim*T_COLS_dim + curr_i_1based*T_COLS_dim + (curr_j_1based-1)
                                m_left  = tl.load(M_matrix_ptr  + left_offset)
                                ix_left = tl.load(Ix_matrix_ptr + left_offset)
                                iy_left = tl.load(Iy_matrix_ptr + left_offset)
                                h_left  = tl.maximum(tl.full((), 0, dtype=tl.int32), tl.maximum(m_left, tl.maximum(ix_left, iy_left)))
                            
                            score_iy_from_h_left = h_left - gap_open

                            iy_prev_left = VERY_NEGATIVE_SCORE_CONST
                            if curr_j_1based > 1 and curr_i_1based >= 1:
                                iy_prev_left_offset = pid*Q_COLS_dim*T_COLS_dim + curr_i_1based*T_COLS_dim + (curr_j_1based-1)
                                iy_prev_left = tl.load(Iy_matrix_ptr + iy_prev_left_offset)
                            
                            score_iy_from_iy_left = iy_prev_left - gap_extend
                            
                            if current_iy_val == score_iy_from_h_left or \
                               (current_iy_val == score_iy_from_iy_left and score_iy_from_h_left < score_iy_from_iy_left):
                                tl.store(out_q_aligned_chars_ptr + pid * MAX_ALIGN_LEN + curr_q_char_write_idx, GAP_CHAR_CODE)
                                tl.store(out_t_aligned_chars_ptr + pid * MAX_ALIGN_LEN + curr_t_char_write_idx, ord_char_t)
                                curr_j_1based -= 1
                                moved_this_step = True
                            elif current_iy_val == score_iy_from_iy_left:
                                tl.store(out_q_aligned_chars_ptr + pid * MAX_ALIGN_LEN + curr_q_char_write_idx, GAP_CHAR_CODE)
                                tl.store(out_t_aligned_chars_ptr + pid * MAX_ALIGN_LEN + curr_t_char_write_idx, ord_char_t)
                                curr_j_1based -= 1
                                moved_this_step = True

                        if moved_this_step:
                            if curr_q_char_write_idx >= 0 and curr_t_char_write_idx >=0 : 
                                curr_q_char_write_idx -= 1
                                curr_t_char_write_idx -= 1
                            current_alignment_len += 1
                        else: 
                            keep_tracing_flag = False
            
            if current_alignment_len > 0:
                tl.store(out_q_aligned_len_ptr + pid, current_alignment_len)
                tl.store(out_t_aligned_len_ptr + pid, current_alignment_len) 
                tl.store(out_q_start_0based_ptr + pid, curr_i_1based) 
                tl.store(out_t_start_0based_ptr + pid, curr_j_1based) 
                tl.store(out_q_end_0based_ptr + pid, q_idx_end_1based - 1) 
                tl.store(out_t_end_0based_ptr + pid, t_idx_end_1based - 1) 
