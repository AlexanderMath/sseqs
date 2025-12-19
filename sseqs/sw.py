'''
python sseqs/sw.py -test_sw
compiled 52.912479639053345
✅ seqs=100000  qlen=256        db_len=50       acc=100/100     TCUPs: 1.27±0.191       time: 0.97±0.15ms       103.08M seqs/s
✅ seqs=100000  qlen=256        db_len=65       acc=100/100     TCUPs: 1.71±0.005       time: 0.93±0.00ms       107.74M seqs/s
✅ seqs=100000  qlen=256        db_len=110      acc=100/100     TCUPs: 2.53±0.006       time: 1.08±0.00ms       92.35M seqs/s
✅ seqs=100000  qlen=256        db_len=129      acc=100/100     TCUPs: 2.23±0.004       time: 1.44±0.00ms       69.24M seqs/s
✅ seqs=100000  qlen=256        db_len=250      acc=100/100     TCUPs: 3.28±0.004       time: 1.93±0.00ms       51.80M seqs/s
✅ seqs=100000  qlen=327        db_len=138      acc=100/100     TCUPs: 2.44±0.001       time: 1.81±0.00ms       55.25M seqs/s
✅ seqs=100000  qlen=77 db_len=67       acc=100/100     TCUPs: 1.44±0.002       time: 0.34±0.00ms       292.82M seqs/s
✅ seqs=100000  qlen=128        db_len=67       acc=100/100     TCUPs: 1.61±0.010       time: 0.51±0.00ms       196.89M seqs/s
✅ seqs=100000  qlen=1024       db_len=1024     acc=100/100     TCUPs: 4.65±0.002       time: 22.47±0.01ms      4.45M seqs/s
✅ seqs=100000  qlen=128        db_len=128      acc=100/100     TCUPs: 2.73±0.003       time: 0.59±0.00ms       170.56M seqs/s
✅ seqs=100000  qlen=1024       db_len=256      acc=100/100     TCUPs: 3.70±0.000       time: 7.00±0.00ms       14.28M seqs/s
✅ seqs=100000  qlen=256        db_len=1024     acc=100/100     TCUPs: 4.22±0.005       time: 6.20±0.01ms       16.13M seqs/s
✅ seqs=100000  qlen=137        db_len=512      acc=100/100     TCUPs: 3.58±0.011       time: 1.95±0.01ms       51.32M seqs/s
✅ seqs=100000  qlen=777        db_len=512      acc=100/100     TCUPs: 4.07±0.000       time: 9.71±0.00ms       10.30M seqs/s
'''
import torch as th 
import triton
import triton.language as tl
import time 
from icecream import ic
import numpy as np 
# wrapper for cuda code, allow using optimized CUDA kernel easily from python. 
from torch.utils.cpp_extension import load
import os 

# Auto-detect GPU and compile only for that architecture
if th.cuda.is_available():
    major, minor = th.cuda.get_device_capability()
    os.environ["TORCH_CUDA_ARCH_LIST"] = f"{major}.{minor}"

package_dir = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()
#sseqs_sw_ext = load(name="sseqs_sw_ext", sources=[f"{package_dir}/sw_bind.cpp", f"{package_dir}/sw.cu"], extra_cuda_cflags=["-O3", "--use_fast_math"],  extra_cflags=["-O3"], verbose=True)
#sseqs_sw_ext = load(name="sseqs_sw_ext", sources=[f"{package_dir}/sw_bind.cpp", f"{package_dir}/sw.cu"], verbose=True)# , extra_cuda_cflags=["--use_fast_math"],  extra_cflags=["-O3"])
sseqs_sw_ext = load("sw",   sources=[f"{package_dir}/sw_bind.cpp",     f"{package_dir}/sw.cu"],     extra_cuda_cflags=["-O3", "--use_fast_math"], extra_cflags=["-O3"], verbose=True)
print('compiled', time.time()-t0)
def sw(query: str, targets: list[str], db_len: int, gap_open=11, gap_extend=1, benchmark=False):
    # Pad targets to multiple of 16 for efficient CUDA processing
    SEQS_PER_BLOCK = 16
    original_len = len(targets)
    remainder = len(targets) % SEQS_PER_BLOCK
    if remainder != 0:
        pad_count = SEQS_PER_BLOCK - remainder
        targets = targets + [targets[-1]] * pad_count
    
    _q_tensor = th.tensor([AA_MAP[aa] for aa in query if aa in AA_MAP], dtype=th.uint8, device="cuda")
    # ASCII concatenation - single string join then one numpy call
    joined = '@'.join(targets) + '@'
    ascii_np = np.frombuffer(joined.encode('latin1'), dtype=np.uint8).copy()  # .copy() makes it writable
    ascii = th.from_numpy(ascii_np).cuda()
    
    delimiter = 64
    starts = th.nonzero(ascii == delimiter, as_tuple=False).flatten()
    starts = th.cat([th.tensor([-1], dtype=th.int32, device="cuda"), starts])
    good_idx = th.arange(1, len(targets)+1, dtype=th.int32, device="cuda")
    th.cuda.synchronize()
    start_event = th.cuda.Event(enable_timing=True)
    end_event = th.cuda.Event(enable_timing=True)
    times = []
    for _ in range(3 if benchmark else 1):
        start_event.record()
        r = sseqs_sw_ext.sw_cuda_affine(
            _q_tensor,
            good_idx,
            ascii,
            starts.to(th.int32),
            db_len,
            gap_open=gap_open,
            gap_extend=gap_extend
        )
        end_event.record()
        th.cuda.synchronize()
        t = start_event.elapsed_time(end_event) 
        times.append(t)
   
    # Return only the original number of results (excluding padding)
    return r[:original_len], times


def sw_linear(query: str, targets: list[str], db_len: int, gap_penalty=1, benchmark=False):
    """Smith-Waterman with linear gap penalties (for screening/prefilter)"""
    SEQS_PER_BLOCK = 16
    original_len = len(targets)
    remainder = len(targets) % SEQS_PER_BLOCK
    if remainder != 0:
        pad_count = SEQS_PER_BLOCK - remainder
        targets = targets + [targets[-1]] * pad_count
    
    _q_tensor = th.tensor([AA_MAP[aa] for aa in query if aa in AA_MAP], dtype=th.uint8, device="cuda")
    joined = '@'.join(targets) + '@'
    ascii_np = np.frombuffer(joined.encode('latin1'), dtype=np.uint8).copy()
    ascii = th.from_numpy(ascii_np).cuda()
    
    delimiter = 64
    starts = th.nonzero(ascii == delimiter, as_tuple=False).flatten()
    starts = th.cat([th.tensor([-1], dtype=th.int32, device="cuda"), starts])
    good_idx = th.arange(1, len(targets)+1, dtype=th.int32, device="cuda")
    th.cuda.synchronize()
    start_event = th.cuda.Event(enable_timing=True)
    end_event = th.cuda.Event(enable_timing=True)
    times = []
    for _ in range(3 if benchmark else 1):
        start_event.record()
        r = sseqs_sw_ext.sw_cuda_linear(
            _q_tensor,
            good_idx,
            ascii,
            starts.to(th.int32),
            db_len,
            gap_penalty=gap_penalty
        )
        end_event.record()
        th.cuda.synchronize()
        t = start_event.elapsed_time(end_event) 
        times.append(t)
   
    return r[:original_len], times


def sw_uint8(query: str, targets: list[str], gap_open=11, gap_extend=1, benchmark=False):
    """Smith-Waterman with uint8 SIMD (for screening, may overflow on high scores)"""
    # Pad targets to multiple of SEQS_PER_BLOCK for int8 kernel
    # (256 threads / 32) * (32 / 8) * 4 = 128 sequences per block (for db_len=128)
    SEQS_PER_BLOCK = 128
    original_len = len(targets)
    remainder = len(targets) % SEQS_PER_BLOCK
    if remainder != 0:
        pad_count = SEQS_PER_BLOCK - remainder
        targets = targets + [targets[-1]] * pad_count
    
    _q_tensor = th.tensor([AA_MAP[aa] for aa in query if aa in AA_MAP], dtype=th.uint8, device="cuda")
    # ASCII concatenation - single string join then one numpy call
    joined = '@'.join(targets) + '@'
    ascii_np = np.frombuffer(joined.encode('latin1'), dtype=np.uint8).copy()  # .copy() makes it writable
    ascii = th.from_numpy(ascii_np).cuda()
    
    delimiter = 64
    starts = th.nonzero(ascii == delimiter, as_tuple=False).flatten()
    starts = th.cat([th.tensor([-1], dtype=th.int32, device="cuda"), starts])
    good_idx = th.arange(1, len(targets)+1, dtype=th.int32, device="cuda")
    th.cuda.synchronize()
    start_event = th.cuda.Event(enable_timing=True)
    end_event = th.cuda.Event(enable_timing=True)
    times = []
    for _ in range(5 if benchmark else 1):
        start_event.record()
        r = sseqs_sw_ext.sw_cuda_affine_uint8(
            _q_tensor,
            good_idx,
            ascii,
            starts.to(th.int32),
            gap_open=gap_open,
            gap_extend=gap_extend
        )
        end_event.record()
        th.cuda.synchronize()
        t = start_event.elapsed_time(end_event) 
        times.append(t)
  
    # Return only the original number of results (excluding padding)
    return r[:original_len], times


def sw_antidiag(query: str, targets: list[str], gap_open=11, gap_extend=1, benchmark=False):
    """ Smith-Waterman with affine gap penalties using anti-diagonal parallelization.
    
    Args:
        query: Query sequence string
        targets: List of target sequence strings
        gap_open: Gap opening penalty (default: 11)
        gap_extend: Gap extension penalty (default: 1)
        benchmark: If True, compute and print TCUPs (default: False)
    
    Returns:
        torch.Tensor of scores for each target sequence
    """
    _q_tensor = th.tensor([AA_MAP[aa] for aa in query if aa in AA_MAP], dtype=th.uint8, device="cuda")
    ascii = th.hstack([th.tensor([ord(c) for c in target+'@'], dtype=th.uint8, device="cuda") for target in targets])
    delimiter = 64
    starts = th.nonzero(ascii == delimiter, as_tuple=False).flatten()
    starts = th.cat([th.tensor([-1], dtype=th.int32, device="cuda"), starts])
    good_idx = th.arange(1, len(targets)+1, dtype=th.int32, device="cuda")
    th.cuda.synchronize()
    for _ in range(2):
        t0 = time.time()
        r = sseqs_sw_ext.sw_cuda_affine_antidiag(
            _q_tensor,
            good_idx,
            ascii,
            starts.to(th.int32),
            gap_open=gap_open,
            gap_extend=gap_extend
        )
        th.cuda.synchronize()
        t = time.time() - t0
        #print(t)
        if benchmark:
            total_cell_updates = len(query) * sum(len(target) for target in targets)
            tcups = total_cell_updates / t / 1e12
            total_tokens = sum(len(target) for target in targets)
            print(f"sw_antidiag TCUPs: {tcups:.2f}, time: {t:.4f}s, sequences: {len(targets)}, total_tokens: {total_tokens}")
    return r, t


def sw_profile(pssm: th.Tensor, targets: list[str], gap_open=11, gap_extend=1, benchmark=False):
    """Smith-Waterman with PSSM (profile) scoring instead of BLOSUM62.
    
    Args:
        pssm: Position-specific scoring matrix, shape (query_len, 20), torch.Tensor (int8, CUDA)
        targets: List of target sequence strings
        gap_open: Gap opening penalty (default: 11)
        gap_extend: Gap extension penalty (default: 1)
        benchmark: If True, run multiple times for timing (default: False)
    
    Returns:
        torch.Tensor of scores for each target sequence, list of times
    """
    import numpy as np
    
    # Ensure PSSM is int8 on CUDA
    if not pssm.is_cuda:
        pssm = pssm.cuda()
    if pssm.dtype != th.int8:
        pssm = pssm.to(th.int8)
    
    query_len = pssm.shape[0]
    
    # Pad targets to multiple of SEQS_PER_BLOCK
    SEQS_PER_BLOCK = 64  # Same as affine kernel
    original_len = len(targets)
    remainder = len(targets) % SEQS_PER_BLOCK
    if remainder != 0:
        pad_count = SEQS_PER_BLOCK - remainder
        targets = targets + [targets[-1]] * pad_count
    
    # Compute max target length for dispatch
    target_len = max(len(t) for t in targets)
    
    # Build concatenated ASCII sequences with delimiter
    joined = '@'.join(targets) + '@'
    ascii_np = np.frombuffer(joined.encode('latin1'), dtype=np.uint8).copy()
    ascii = th.from_numpy(ascii_np).cuda()
    
    delimiter = 64  # '@' 
    starts = th.nonzero(ascii == delimiter, as_tuple=False).flatten()
    starts = th.cat([th.tensor([-1], dtype=th.int32, device="cuda"), starts])
    good_idx = th.arange(1, len(targets)+1, dtype=th.int32, device="cuda")
    
    th.cuda.synchronize()
    start_event = th.cuda.Event(enable_timing=True)
    end_event = th.cuda.Event(enable_timing=True)
    times = []
    
    for _ in range(3 if benchmark else 1):
        start_event.record()
        r = sseqs_sw_ext.sw_cuda_profile(
            pssm,
            good_idx,
            ascii,
            starts.to(th.int32),
            target_len,
            gap_open=gap_open,
            gap_extend=gap_extend
        )
        end_event.record()
        th.cuda.synchronize()
        t = start_event.elapsed_time(end_event) / 1000.0  # Convert ms to seconds
        times.append(t)
    
    return r[:original_len], times


def sw_affine_backtrack_cuda(query: str, targets: list[str], gap_open=11, gap_extend=1):
    """CUDA implementation with checkpointing for memory-efficient backtracking."""
    _q_tensor = th.tensor([AA_MAP[aa] for aa in query if aa in AA_MAP], dtype=th.uint8, device="cuda")
    ascii = th.hstack([th.tensor([ord(c) for c in target+'@'], dtype=th.uint8, device="cuda") for target in targets])
    delimiter = 64
    starts = th.nonzero(ascii == delimiter, as_tuple=False).flatten()
    starts = th.cat([th.tensor([-1], dtype=th.int32, device="cuda"), starts])
    good_idx = th.arange(1, len(targets)+1, dtype=th.int32, device="cuda")
    
    max_align_len = len(query) + max(len(t) for t in targets)
    import time 
    for _ in range(1):
        t0 = time.time()
        scores, align_lens, end_i, end_j, align_ops = sseqs_sw_ext.sw_cuda_affine_backtrack(
            _q_tensor,
            good_idx,
            ascii,
            starts.to(th.int32),
            max_align_len,
            gap_open,
            gap_extend
        )
        th.cuda.synchronize()
        t= time.time()-t0
        ic(t)
    
    # Reconstruct aligned sequences from operations
    alignments = []
    for i in range(len(targets)):
        alen = align_lens[i].item()
        if alen == 0:
            alignments.append({'q_aligned': '', 't_aligned': '', 'score': scores[i].item()})
            continue
            
        ops = align_ops[i, :alen].cpu().numpy().tobytes().decode('latin-1')
        
        # CUDA DP indexing: When at cell (i, j) and we take a MATCH, we consume target[i] and query[j]
        # The end positions from CUDA are the (i, j) of the last cell in the alignment
        # So the last MATCH consumed target[end_i] and query[end_j]
        # 
        # Working backwards through the operations:
        # - Each MATCH: consume both target and query (i--, j--)
        # - Each INSERT: consume only query (j--)
        # - Each DELETE: consume only target (i--)
        #
        # After we've gone through all operations backward, we're at the cell BEFORE the first operation
        # So we need to add 1 to get the index of the first character consumed
        
        # Calculate the position before the first operation
        q_before, t_before = end_j[i].item(), end_i[i].item()
        for op in ops:
            if op == 'M':
                q_before -= 1
                t_before -= 1
            elif op == 'I':
                q_before -= 1
            elif op == 'D':
                t_before -= 1
        
        # The start position is one after the "before" position
        q_start, t_start = q_before + 1, t_before + 1
        
        # DEBUG: Print for first sequence
        if i == 0 and False:  # Disable debug output
            print(f"\n=== DEBUG Seq {i} ===")
            print(f"Target: {targets[i]}")
            print(f"Query:  {query}")
            print(f"Ops: {ops}")
            print(f"End pos: t={end_i[i].item()} q={end_j[i].item()}")
            print(f"Start pos: t={t_start} q={q_start}")
            print(f"Target substring [{t_start}:{end_i[i].item()}]: {targets[i][t_start:end_i[i].item()]}")
            print(f"Query substring [{q_start}:{end_j[i].item()}]: {query[q_start:end_j[i].item()]}")
            print(f"\nTrying different interpretations:")
            print(f"  [{t_start}:{end_i[i].item()+1}] vs [{q_start}:{end_j[i].item()+1}]:")
            print(f"    Target: {targets[i][t_start:end_i[i].item()+1]}")
            print(f"    Query:  {query[q_start:end_j[i].item()+1]}")
            print(f"  [{t_start+1}:{end_i[i].item()+1}] vs [{q_start+1}:{end_j[i].item()+1}]:")
            print(f"    Target: {targets[i][t_start+1:end_i[i].item()+1]}")
            print(f"    Query:  {query[q_start+1:end_j[i].item()+1]}")
        
        # Now build alignment going FORWARD from start position
        q_aln, t_aln = [], []
        q_pos, t_pos = q_start, t_start
        
        for op in ops:
            if op == 'M':  # Match/mismatch
                if q_pos >= 0 and t_pos >= 0 and q_pos < len(query) and t_pos < len(targets[i]):
                    q_aln.append(query[q_pos])
                    t_aln.append(targets[i][t_pos])
                q_pos += 1
                t_pos += 1
            elif op == 'I':  # Insertion in query (gap in target)
                if q_pos >= 0 and q_pos < len(query):
                    q_aln.append(query[q_pos])
                    t_aln.append('-')
                q_pos += 1
            elif op == 'D':  # Deletion in query (gap in query)
                if t_pos >= 0 and t_pos < len(targets[i]):
                    q_aln.append('-')
                    t_aln.append(targets[i][t_pos])
                t_pos += 1
        
        if i == 0:
            print(f"Extracted q_aln: {repr(''.join(q_aln))}")
            print(f"Extracted t_aln: {repr(''.join(t_aln))}")
        
        alignments.append({
            'q_aligned': ''.join(q_aln),
            't_aligned': ''.join(t_aln),
            'score': scores[i].item()
        })
    
    return scores, alignments, t

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

eval_file = f'eval/runs_and_poses/raw_us7/raw_us7/7fak__1__1.A__1.B.pkl'
eval_file = f"eval/runs_and_poses/low_homology_raw_all/553a03b7ac0b44b7a6293d3c90bfae1df25d04291fe784db5b38f87afcb88344.pkl"

def test_sw(profile):
    import pickle 
    import numpy as np 
    import parasail
    import torch
    np.random.seed(42)
    # if we could train a int4 NN we could get 60 TCUPs; this is faster than PCIE i guess? 
    
    # wget foldify.org/test_case.pkl #  80mb of hits
    q, seqs, scores = pickle.load(open(eval_file, 'rb'))
    seqs = [s.decode('ascii') for s in seqs]
    q = q * 3
    seqs = sorted(seqs, key=len)[::-1]
    seqs = [s*10 for s in seqs[:100_000]] 
    seqs = [s.replace('@','').replace('B', '').replace('X', '').replace('-','').replace('*','') for s in seqs if len(s)>1024]
   
    for num, (target_len, query_len) in enumerate([
        [50,       256],
        [65,       256],
        [110,      256],
        [129,      256],
        [250,      256],
        [138,       327],
        [67,       77],
        [67,       128],
        [1024,     1024],
        #[5,         7],
        #[129,       128],
        [128,       128],
        [256,       1024],
        [1024,      256],
        [512,      137],
        [512,      777],
       ]):
        # Prepare sequences
        targets = [s[:target_len+np.random.randint(-5, 0)] for s in seqs]
        query = q[:query_len]
        
        # Benchmark CUDA on all sequences
        scores_cuda, times = sw(query, targets, target_len, benchmark=True, gap_open=11, gap_extend=1)
        
        # Compute TCUPs stats (skip first cold run)
        times_arr = np.array(times[1:]) / 1000.0  # Convert ms to seconds for TCUPs
        median_t = np.median(times_arr)
        std_t = np.std(times_arr)
        total_cell_updates = len(query) * sum(len(target) for target in targets)
        tcups_arr = total_cell_updates / times_arr / 1e12
        median_tcups = np.median(tcups_arr)
        std_tcups = np.std(tcups_arr)
        
        # Validate on random sample
        n_validate = 100
        validate_indices = np.random.choice(len(targets), size=min(n_validate, len(targets)), replace=False)
        scores_parasail = []
        for i in validate_indices:
            result = parasail.sw_trace(query, targets[i], 11, 1, parasail.blosum62)
            scores_parasail.append(result.score)
        scores_parasail = torch.tensor(scores_parasail, device='cuda')
        
        # Check accuracy (allow tolerance for FP16 precision in large alignments)
        scores_cuda_sampled = scores_cuda[validate_indices]
        diff = (scores_cuda_sampled - scores_parasail).abs()
        rel_error = diff / scores_parasail.float()
        exact_matches = (rel_error < 0.03).sum().item()
        acc_emoji = "✅" if exact_matches == n_validate else "❌"
        
        # If not 100% accurate, report errors
        if exact_matches != n_validate:
            wrong_indices = (rel_error >= 0.03).nonzero(as_tuple=True)[0]
            print(f"  Errors found in {len(wrong_indices)} sequences:")
            for idx in wrong_indices[:10]:  # Show first 10 errors
                i = idx.item()
                cuda_score = scores_cuda_sampled[i].item()
                parasail_score = scores_parasail[i].item()
                error = rel_error[i].item()
                print(f"    idx={validate_indices[i]}: CUDA={cuda_score:.1f} Parasail={parasail_score:.1f} rel_err={error:.3f}")
            if len(wrong_indices) > 10:
                print(f"    ... and {len(wrong_indices) - 10} more errors")
        
        # Compute sequences per second
        seqs_per_sec = len(targets) / median_t / 1e6  # In millions
        
        # Single line output (times back to ms for display)
        print(f"{acc_emoji} seqs={len(targets)}\tqlen={query_len}\tdb_len={target_len}\tacc={exact_matches}/{n_validate}\tTCUPs: {median_tcups:.2f}±{std_tcups:.3f}\ttime: {median_t*1000:.2f}±{std_t*1000:.2f}ms\t{seqs_per_sec:.2f}M seqs/s")
        if profile: exit()


def test_sw_backtrack(profile):
    """Test SW with backtracking"""
    import pickle 
    import numpy as np 
    import parasail
    import torch 
    np.random.seed(42)
    # 5TCUPs = 4m seqs/s.  i just want ~2m for backtrack.
    
    print("Testing SW Backtrack (sw_affine_backtrack_cuda)")
    print("=" * 80)
    
    # Load test data
    q, seqs, scores = pickle.load(open(eval_file, 'rb'))
    seqs = [seq.decode('ascii') for seq in seqs]
    seqs = [s.replace('@','').replace('B', '').replace('X', '').replace('-','').replace('*','') for s in seqs if len(s)>50]
    seqs = [s[:1024] for s in seqs]
    
    # Sort by length (longest first)
    seqs = sorted(seqs, key=len, reverse=True)
    seqs = seqs[:10]  # Test on 10 sequences only
    query = q[:256]
    from icecream import ic 
    ic(seqs)
    ic(query)
    
    print(f"Query length: {len(query)}")
    print(f"Testing on {len(seqs)} sequences")
    print(f"Target lengths: {[len(s) for s in seqs[:10]]}")
    
    # First test: Forward pass scores should match parasail
    print("\n=== Test 1: Forward Pass Scores ===")
    print("Comparing backtrack kernel scores vs parasail...")
    
    # Run CUDA backtrack kernel
    scores_backtrack, alignments, t = sw_affine_backtrack_cuda(query, seqs)
    print(t)
    
    # Get parasail reference scores
    import parasail
    scores_parasail = []
    alignments_parasail = []
    for seq in seqs:
        result = parasail.sw_trace(query, seq, 11, 1, parasail.blosum62)
        scores_parasail.append(result.score)
        
        # Get the actual alignment from parasail for comparison
        if len(alignments_parasail) < 1:  # Just first one for debug
            cigar = result.cigar
            alignments_parasail.append({
                'score': result.score,
                'end_ref': result.end_ref,
                'end_query': result.end_query,
                'len_ref': result.len_ref,
                'len_query': result.len_query,
                'cigar': cigar.decode if hasattr(cigar, 'decode') else str(cigar)
            })
    scores_parasail = torch.tensor(scores_parasail, device='cuda')
    
    print(f"Parasail scores:     {scores_parasail.cpu().numpy()[:10]}")
    print(f"Backtrack SW scores: {scores_backtrack.cpu().numpy()[:10]}")
    
    # DEBUG: Show parasail's alignment for first sequence
    if alignments_parasail:
        p = alignments_parasail[0]
        print(f"\nParasail alignment for seq 0:")
        print(f"  Score: {p['score']}")
        print(f"  end_ref: {p['end_ref']}, end_query: {p['end_query']}")
        print(f"  len_ref: {p['len_ref']}, len_query: {p['len_query']}")
        ref_begin = p['end_ref'] - p['len_ref'] + 1
        query_begin = p['end_query'] - p['len_query'] + 1
        print(f"  Computed: ref_begin={ref_begin}, query_begin={query_begin}")
        print(f"  Target: [{ref_begin}:{p['end_ref']+1}] = {seqs[0][ref_begin:p['end_ref']+1]}")
        print(f"  Query:  [{query_begin}:{p['end_query']+1}] = {query[query_begin:p['end_query']+1]}")
        print(f"  CIGAR: {p['cigar']}")
    
    # Compare scores (allow 2% relative difference)
    score_diff = (scores_backtrack - scores_parasail).abs()
    relative_diff = score_diff / (scores_parasail.abs() + 1e-6)  # Avoid division by zero
    max_rel_diff = relative_diff.max().item()
    max_abs_diff = score_diff.max().item()
    matches = (relative_diff <= 0.02).sum().item()
    
    if matches == len(seqs):
        print(f"✅ Forward pass: ALL scores match parasail within 2%! (max rel diff: {max_rel_diff*100:.2f}%, max abs diff: {max_abs_diff:.1f})")
    else:
        print(f"❌ Forward pass: {matches}/{len(seqs)} scores match parasail (max rel diff: {max_rel_diff*100:.2f}%)")
        for i in range(len(seqs)):
            if relative_diff[i] > 0.02:
                rel_pct = relative_diff[i].item() * 100
                print(f"  Seq {i}: Parasail={scores_parasail[i].item():.0f} Backtrack={scores_backtrack[i].item():.0f} diff={score_diff[i].item():.0f} ({rel_pct:.2f}%)")
        print("\n❌ Forward pass scores don't match - fix the kernel before testing backtracking!")
        return
    
    # Second test: Validate backtracking correctness
    print("\n=== Test 2: Backtracking Correctness ===")
    print("Testing: Does parasail agree our alignment achieves the claimed score?")
    
    for i in range(len(seqs)):
        cuda_score = scores_backtrack[i].item()
        
        # Get pre-aligned sequences from CUDA backtrack
        q_aligned = alignments[i]['q_aligned']
        t_aligned = alignments[i]['t_aligned']
        
        # Remove gaps to get the actual subsequences that were aligned
        q_nogap = q_aligned.replace('-', '')
        t_nogap = t_aligned.replace('-', '')
        
        # Re-score using SW on the gap-free subsequences
        result_rescore = parasail.sw_trace(q_nogap, t_nogap, 11, 1, parasail.blosum62)
        parasail_rescore = result_rescore.score
        
        # Check if parasail agrees this alignment has the claimed score (allow 2% relative difference)
        diff = abs(cuda_score - parasail_rescore)
        rel_diff = diff / (abs(cuda_score) + 1e-6)
        status = "✅" if rel_diff <= 0.02 else "❌"
        
        print(f"{status} Seq {i}: CUDA_score={cuda_score:.0f} Parasail_rescore={parasail_rescore:.0f} diff={diff:.0f} ({rel_diff*100:.2f}%) align_len={len(q_aligned)}")
        
        # Show details for failures
        if rel_diff > 0.02:
            print(f"  Query nogap:  {q_nogap[:30]}...")
            print(f"  Target nogap: {t_nogap[:30]}...")
            print(f"  Query aligned:  {q_aligned[:30]}...")
            print(f"  Target aligned: {t_aligned[:30]}...")
    
    print("\n✅ Backtrack test complete!")


def test_sw_linear(profile):
    """Test linear gap penalty kernel"""
    import pickle 
    import numpy as np 
    import parasail
    import torch
    np.random.seed(42)
    
    print("Testing Linear Gap Penalty Kernel (sw_kernel_linear)")
    print("=" * 80)
    
    # Load test data
    q, seqs, scores = pickle.load(open(eval_file, 'rb'))
    q = q * 2
    seqs = sorted(seqs, key=len)[::-1]
    seqs = [s*5 for s in seqs[:100_000]] 
    seqs = [s.replace('@','').replace('B', '').replace('X', '').replace('-','').replace('*','') for s in seqs if len(s)>128]
   
    for num, (target_len, query_len) in enumerate([
        [128,      128],
        [256,      256],
        [512,      512],
        [1024,     1024],
       ]):
        # Prepare sequences
        targets = [s[:target_len+np.random.randint(-5, 0)] for s in seqs]
        query = q[:query_len]
        
        # Benchmark CUDA linear kernel
        scores_cuda_linear, times = sw_linear(query, targets, target_len, benchmark=True, gap_penalty=1)
        
        # Compute TCUPs stats
        times_arr = np.array(times[1:]) / 1000.0
        median_t = np.median(times_arr)
        std_t = np.std(times_arr)
        total_cell_updates = len(query) * sum(len(target) for target in targets)
        tcups_arr = total_cell_updates / times_arr / 1e12
        median_tcups = np.median(tcups_arr)
        std_tcups = np.std(tcups_arr)
        
        # Validate against parasail with linear gaps
        # Note: parasail doesn't have pure linear gaps, so we use affine with gap_open = gap_extend
        n_validate = min(100, len(targets))
        validate_indices = np.random.choice(len(targets), size=n_validate, replace=False)
        scores_parasail = []
        for i in validate_indices:
            # Simulate linear gaps: set gap_open = gap_extend = 1
            result = parasail.sw_trace(query, targets[i], 1, 1, parasail.blosum62)
            scores_parasail.append(result.score)
        scores_parasail = torch.tensor(scores_parasail, device='cuda')
        
        # Check accuracy
        scores_cuda_sampled = scores_cuda_linear[validate_indices]
        diff = (scores_cuda_sampled - scores_parasail).abs()
        rel_error = diff / (scores_parasail.float() + 1e-6)
        exact_matches = (rel_error < 0.03).sum().item()
        acc_emoji = "✅" if exact_matches == n_validate else "❌"
        
        # Report any errors
        if exact_matches != n_validate:
            wrong_indices = (rel_error >= 0.03).nonzero(as_tuple=True)[0]
            print(f"  Errors found in {len(wrong_indices)} sequences:")
            for idx in wrong_indices[:5]:
                i = idx.item()
                cuda_score = scores_cuda_sampled[i].item()
                parasail_score = scores_parasail[i].item()
                error = rel_error[i].item()
                print(f"    idx={validate_indices[i]}: CUDA={cuda_score:.1f} Parasail={parasail_score:.1f} rel_err={error:.3f}")
        
        seqs_per_sec = len(targets) / median_t / 1e6
        
        print(f"{acc_emoji} LINEAR seqs={len(targets)}\tqlen={query_len}\tdb_len={target_len}\tacc={exact_matches}/{n_validate}\tTCUPs: {median_tcups:.2f}±{std_tcups:.3f}\ttime: {median_t*1000:.2f}±{std_t*1000:.2f}ms\t{seqs_per_sec:.2f}M seqs/s")
        if profile: exit()

def sw_profile_reference(pssm, target, gap_open=11, gap_extend=1):
    """Pure Python reference implementation for profile SW (slow but correct).
    
    Args:
        pssm: numpy array (query_len, 20) - position-specific scoring matrix
        target: string - target sequence
        gap_open: gap opening penalty
        gap_extend: gap extension penalty
    
    Returns:
        Maximum SW score (int)
    """
    AA_ORDER = 'ARNDCQEGHILKMFPSTWYV'
    aa_to_idx = {aa: i for i, aa in enumerate(AA_ORDER)}
    
    L = len(pssm)  # query length
    T = len(target)
    
    # DP matrices
    H = np.zeros((L + 1, T + 1), dtype=np.float32)
    E = np.zeros((L + 1, T + 1), dtype=np.float32)  # horizontal gap (insertion in query)
    F = np.zeros((L + 1, T + 1), dtype=np.float32)  # vertical gap (deletion in query)
    
    max_score = 0
    for i in range(1, L + 1):
        for j in range(1, T + 1):
            # Get target amino acid index
            t_char = target[j - 1]
            t_idx = aa_to_idx.get(t_char, 0)  # default to 0 (Ala) for unknown
            
            # Match/mismatch score from PSSM
            match = H[i-1, j-1] + pssm[i-1, t_idx]
            
            # Gap scores (affine)
            E[i, j] = max(E[i, j-1] - gap_extend, H[i, j-1] - gap_open)
            F[i, j] = max(F[i-1, j] - gap_extend, H[i-1, j] - gap_open)
            
            # Cell score
            H[i, j] = max(0, match, E[i, j], F[i, j])
            max_score = max(max_score, H[i, j])
    
    return int(max_score)


def test_sw_profile(profile):
    """Test Profile-based Smith-Waterman (PSSM scoring)"""
    import pickle 
    import numpy as np 
    import torch 
    np.random.seed(42)
    
    print("Testing Profile SW (sw_profile)")
    print("=" * 80)
    
    # Load test data
    q, seqs, scores = pickle.load(open(eval_file, 'rb'))
    seqs = [seq.decode('ascii') for seq in seqs]
    seqs = [s.replace('@','').replace('B', '').replace('X', '').replace('-','').replace('*','') for s in seqs if len(s) > 50]
    
    # Use smaller subset for testing
    query = q[:128]  # shorter for faster Python reference
    query_len = len(query)
    
    # Sort by length and filter to reasonable size
    seqs = sorted(seqs, key=len, reverse=True)
    target_len = 128  # shorter for faster Python reference
    targets = [s[:target_len] for s in seqs if len(s) >= target_len][:10000]
    
    print(f"Query length: {query_len}")
    print(f"Testing on {len(targets)} sequences")
    print(f"Target length: {target_len}")
    
    AA_ORDER = 'ARNDCQEGHILKMFPSTWYV'
    aa_to_idx = {aa: i for i, aa in enumerate(AA_ORDER)}
    
    blosum62 = np.array([
        [ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0],  # A
        [-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3],  # R
        [-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3],  # N
        [-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3],  # D
        [ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],  # C
        [-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2],  # Q
        [-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2],  # E
        [ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3],  # G
        [-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3],  # H
        [-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3],  # I
        [-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1],  # L
        [-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2],  # K
        [-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1],  # M
        [-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1],  # F
        [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2],  # P
        [ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2],  # S
        [ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0],  # T
        [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3],  # W
        [-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1],  # Y
        [ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4],  # V
    ], dtype=np.int8)
    
    # =====================================================
    # Test 1: BLOSUM62-derived PSSM vs regular affine SW
    # =====================================================
    print("\n=== Test 1: BLOSUM62-derived PSSM vs Affine SW ===")
    
    # Build PSSM from query using BLOSUM62
    pssm_blosum = np.zeros((query_len, 20), dtype=np.int8)
    for i, q_char in enumerate(query):
        if q_char in aa_to_idx:
            q_idx = aa_to_idx[q_char]
            pssm_blosum[i, :] = blosum62[q_idx, :]
        else:
            pssm_blosum[i, :] = 0
    
    print(f"BLOSUM-derived PSSM shape: {pssm_blosum.shape}")
    
    pssm_blosum_tensor = torch.from_numpy(pssm_blosum).cuda()
    scores_profile_blosum, _ = sw_profile(pssm_blosum_tensor, targets)
    scores_affine, _ = sw(query, targets, db_len=target_len)
    
    diff_blosum = (scores_profile_blosum - scores_affine).abs()
    max_diff_blosum = diff_blosum.max().item()
    matches_blosum = (diff_blosum == 0).sum().item()
    
    print(f"Exact matches: {matches_blosum}/{len(targets)}")
    print(f"Max absolute difference: {max_diff_blosum}")
    
    # Print sample scores to verify they're not garbage
    print(f"Sample scores (first 10): {scores_profile_blosum[:10].cpu().numpy()}")
    
    if max_diff_blosum == 0:
        print(f"✅ Profile SW matches Affine SW with BLOSUM-derived PSSM!")
    else:
        print(f"❌ Differences found:")
        wrong_idx = (diff_blosum > 0).nonzero(as_tuple=True)[0][:5]
        for idx in wrong_idx:
            i = idx.item()
            print(f"  idx={i}: Profile={scores_profile_blosum[i].item()} Affine={scores_affine[i].item()}")
    
    # =====================================================
    # Test 2: Random PSSM vs Python reference (cached)
    # =====================================================
    print("\n=== Test 2: Random PSSM vs Python Reference ===")
    
    # Generate random PSSM
    np.random.seed(42)
    pssm_random = np.random.randint(-4, 12, size=(query_len, 20)).astype(np.int8)
    print(f"Random PSSM shape: {pssm_random.shape}, range: [{pssm_random.min()}, {pssm_random.max()}]")
    
    n_validate = 100
    validate_targets = targets[:n_validate]
    
    # Check for cached reference scores
    cache_file = f"{package_dir}/profile_ref.pkl"
    cache_valid = False
    
    if os.path.exists(cache_file):
        try:
            cached = pickle.load(open(cache_file, 'rb'))
            # Verify cache matches current test setup
            if (cached['query_len'] == query_len and 
                cached['target_len'] == target_len and 
                cached['n_validate'] == n_validate and
                np.array_equal(cached['pssm'], pssm_random) and
                cached['targets'] == validate_targets):
                scores_ref = cached['scores']
                cache_valid = True
                print(f"Loaded cached reference scores from {cache_file}")
        except Exception as e:
            print(f"Cache invalid: {e}")
    
    if not cache_valid:
        print(f"Computing Python reference scores for {n_validate} sequences...")
        scores_ref = []
        for i, target in enumerate(validate_targets):
            score = sw_profile_reference(pssm_random, target)
            scores_ref.append(score)
            if (i + 1) % 20 == 0:
                print(f"  {i+1}/{n_validate}...")
        scores_ref = np.array(scores_ref)
        
        # Save to cache
        cache_data = {
            'query_len': query_len,
            'target_len': target_len,
            'n_validate': n_validate,
            'pssm': pssm_random,
            'targets': validate_targets,
            'scores': scores_ref
        }
        pickle.dump(cache_data, open(cache_file, 'wb'))
        print(f"Saved reference scores to {cache_file}")
    
    # Run CUDA
    pssm_random_tensor = torch.from_numpy(pssm_random).cuda()
    scores_cuda, _ = sw_profile(pssm_random_tensor, validate_targets)
    scores_cuda_np = scores_cuda.cpu().numpy()
    
    # Compare
    diff_ref = np.abs(scores_cuda_np - scores_ref)
    max_diff_ref = diff_ref.max()
    exact_matches_ref = (diff_ref == 0).sum()
    
    print(f"Exact matches: {exact_matches_ref}/{n_validate}")
    print(f"Max absolute difference: {max_diff_ref}")
    
    # Print sample scores to verify they're not garbage
    print(f"Sample CUDA scores (first 10): {scores_cuda_np[:10]}")
    print(f"Sample Ref  scores (first 10): {scores_ref[:10]}")
    
    if exact_matches_ref == n_validate:
        print(f"✅ CUDA matches Python reference perfectly!")
    else:
        print(f"❌ Differences found:")
        wrong_idx = np.where(diff_ref > 0)[0][:10]
        for idx in wrong_idx:
            print(f"  idx={idx}: CUDA={scores_cuda_np[idx]} Ref={scores_ref[idx]} diff={diff_ref[idx]}")
    
    # =====================================================
    # Test 3: Performance benchmark
    # =====================================================
    print("\n=== Test 3: Performance ===")
    scores_cuda, t = sw_profile(pssm_random_tensor, targets, benchmark=True)
    
    median_t = np.median(t)
    tcups = (query_len * target_len * len(targets)) / median_t / 1e12
    seqs_per_sec = len(targets) / median_t / 1e6
    
    print(f"Profile SW: TCUPs: {tcups:.2f}, {seqs_per_sec:.2f}M seqs/s, {median_t*1000:.2f}ms")
    
    if profile:
        exit()



if __name__ == '__main__':
    import argparse 
    parser = argparse.ArgumentParser(description="Run sequence alignment using SW algorithm.")
    parser.add_argument("-profile", action="store_true", help="Profile the execution time.")
    parser.add_argument("-n",  default=100_000, type=int, help="Number of sequences to test. ")
    parser.add_argument("-test_sw",  action='store_true', help="Test SW used for search. ")
    parser.add_argument("-test_swl",  action='store_true', help="Test SW Linear (screening kernel). ")
    parser.add_argument("-test_swb",  action='store_true', help="Test SW affine backtrack (alignment kernel). ")
    parser.add_argument("-test_sw_uint8",  action='store_true', help="Test SW linear uint8 (screening kernel). ")
    parser.add_argument("-test_swp",  action='store_true', help="Test SW profile (screening kernel). ")
    args = parser.parse_args()

    if args.test_sw: test_sw(args.profile); exit()
    if args.test_swl: test_sw_linear(args.profile); exit()
    if args.test_swb: test_sw_backtrack(args.profile); exit()
    if args.test_sw_uint8: test_sw_uint8(args.profile); exit()
    if args.test_swp: test_sw_profile(args.profile); exit()

    