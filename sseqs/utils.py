import torch
import numpy as np

eval_file = 'eval/runs_and_poses/mmseqs_uniref_plus_bfd/553a03b7ac0b44b7a6293d3c90bfae1df25d04291fe784db5b38f87afcb88344.a3m'

def npz2msa(npz_path):
    """
    Convert Boltz NPZ MSA format to list of sequences.
    
    Args:
        npz_path: Path to the NPZ file containing MSA data
    
    Returns:
        list of sequences (strings with gaps '-')
    """
    # Boltz alphabet: gap + 20 amino acids + unknown (1-indexed in NPZ)
    alphabet = "-ARNDCQEGHILKMFPSTWYVX"
    
    # Load NPZ file
    x = np.load(npz_path, allow_pickle=True)
    meta = x['sequences']
    res = x['residues']
    
    # Extract sequences using the alphabet
    # res_type values are 1-indexed, so we subtract 1 to get 0-indexed alphabet
    seqs = [
        "".join(alphabet[int(r) - 1] for r in res[int(m[2]):int(m[3])]['res_type'])
        for m in meta
    ]
    
    return seqs

def msa_stats_npz(path, num=2048):
    return msa_stats_seqs(npz2msa(path)[:num]) 


def neff(path, threshold=0.62, max_lines=1024):
    """
    Compute N_eff (effective number of sequences). 
    
    Formula: N_eff = Σ_n (1 / num_cluster_members_n)
    
    Where for each sequence n:
    - Identity(i,j) = matches(i,j) / non_gap_positions(i)
    - matches(i,j) = positions where both i and j are non-gap AND equal
    - num_cluster_members_n = number of sequences with Identity >= threshold (including self)
    
    This excludes gaps from the identity calculation 
    
    Args:
        path: Path to MSA file (a3m format) or npz file
        threshold: Sequence identity threshold (default 0.62, typical range 0.62-0.67)
        max_lines: Maximum number of sequences to process
    
    Returns:
        float: Effective number of sequences (N_eff)
    """
    # Load sequences based on file type
    if path.endswith('.npz'):
        sequences = npz2msa(path)
    else:
        lines = open(path, 'r').read().split('\n')
        sequences = lines[2:][1::2]
    
    # Remove lowercase letters (insertions relative to query)
    remove_lower = str.maketrans('', '', 'abcdefghijklmnopqrstuvwxyz')
    seqs_clean = [seq.translate(remove_lower) for seq in sequences if seq.strip()]
    
    if len(seqs_clean) == 0:
        return 0.0

    if len(seqs_clean) == 0:
        return 0.0

    seqs_clean = seqs_clean[:max_lines] 
    
    if len(seqs_clean) == 0:
        return 0.0
    
    # In a3m format (after removing lowercase), all sequences have the same length (aligned)
    seq_len = len(seqs_clean[0])
    M = len(seqs_clean)
    
    # Create array directly - no padding needed
    a3m = np.array([list(seq) for seq in seqs_clean], dtype='U1')
    
    num_sequences, seq_len = a3m.shape
    
    # Use GPU if available for faster computation
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # Convert string array to integer encoding
    unique_chars = np.unique(a3m.ravel())
    char_to_int = {char: i for i, char in enumerate(unique_chars)}
    
    # Encode as integers
    a3m_int = np.vectorize(char_to_int.get)(a3m).astype(np.uint8)
    gap_code = char_to_int['-']
    
    # Convert to torch tensor and move to GPU
    a3m_tensor = torch.from_numpy(a3m_int).to(device).to(torch.uint8)
    not_gap = (a3m_tensor != gap_code)
    
    # Compute sequence lengths (non-gap positions per sequence)
    # Identity is computed as: matches / non_gap_length[i] for sequence i
    seq_lengths = not_gap.sum(dim=1, dtype=torch.float32)  # (num_seqs,)
    
    # Compute pairwise sequence identities (exclude gaps)
    a3m_i = a3m_tensor.unsqueeze(1)  # (num_seqs, 1, seq_len)
    a3m_j = a3m_tensor.unsqueeze(0)  # (1, num_seqs, seq_len)
    
    not_gap_i = not_gap.unsqueeze(1)  # (num_seqs, 1, seq_len)
    not_gap_j = not_gap.unsqueeze(0)  # (1, num_seqs, seq_len)
    
    # Matches: same character AND both non-gap (gaps excluded from identity)
    both_not_gap = not_gap_i & not_gap_j
    matches = (a3m_i == a3m_j) & both_not_gap
    
    # Count matches per pair
    num_matches = matches.sum(dim=2, dtype=torch.float32)  # (num_seqs, num_seqs)
    
    # Compute identities: matches / non_gap_length[i] for each row i
    # This makes the identity matrix asymmetric (Identity[i,j] != Identity[j,i])
    identities = num_matches / seq_lengths.unsqueeze(1).clamp(min=1)  # (num_seqs, num_seqs)
    
    # For each sequence i, count how many sequences (including self) have identity >= threshold
    above_threshold = (identities >= threshold).to(torch.float32)  # (N, N)
    num_cluster_members = above_threshold.sum(dim=1)  # (N,) - includes self
    
    # Compute weights: 1 / num_cluster_members
    weights = 1.0 / num_cluster_members
    
    # N_eff: sum of weights
    neff_value = weights.sum().cpu().item()
    
    # Clean up GPU memory
    del a3m_tensor, not_gap, a3m_i, a3m_j, not_gap_i, not_gap_j, both_not_gap, matches
    torch.cuda.empty_cache()
    
    return float(neff_value)


def msa_stats(path):
    """
    Compute N_eff and coverage for an MSA - FAST GPU version.
    
    Args:
        sequences: list of aligned sequences (strings with gaps '-')
    
    Returns:
        dict with N_eff, coverage, and other stats
    """
    lines = open(path, 'r').read().split('\n')
    sequences = lines[2:][1::2]
    return msa_stats_seqs(sequences)

def msa_stats_seqs(sequences):
    """
    Compute N_eff and coverage for an MSA - FAST GPU version.
    
    Args:
        sequences: list of aligned sequences (strings with gaps '-')
    
    Returns:
        dict with N_eff, coverage, and other stats
    """

    # Remove lowercase letters (insertions relative to query)
    remove_lower = str.maketrans('', '', 'abcdefghijklmnopqrstuvwxyz')
    seqs_clean = [seq.translate(remove_lower) for seq in sequences]
    
    # Find max length and pad all sequences
    max_len = max(len(s) for s in seqs_clean)
    M = len(seqs_clean)
    
    # Create padded array
    a3m = np.full((M, max_len), '-', dtype='U1')
    for i, seq in enumerate(seqs_clean):
        a3m[i, :len(seq)] = list(seq)
    
    num_sequences, seq_len = a3m.shape
    
    # Average length (excluding gaps)
    non_gap_counts = np.sum(a3m != '-', axis=1)
    avg_length = float(np.mean(non_gap_counts))
    
    # Move to GPU for fast pairwise computations
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # Convert string array to integer encoding: '-' = 0, 'A' = 1, etc.
    unique_chars = np.unique(a3m.ravel())
    char_to_int = {char: i for i, char in enumerate(unique_chars)}
    
    # Encode as integers - use uint8 to save memory
    a3m_int = np.vectorize(char_to_int.get)(a3m).astype(np.uint8)
    gap_code = char_to_int['-']
    
    # Convert to torch tensor and move to GPU
    a3m_tensor = torch.from_numpy(a3m_int).to(device).to(torch.uint8)
    not_gap = (a3m_tensor != gap_code)
    
    # Pairwise identity computation on GPU
    a3m_i = a3m_tensor.unsqueeze(1)  # (num_seqs, 1, seq_len)
    a3m_j = a3m_tensor.unsqueeze(0)  # (1, num_seqs, seq_len)
    
    not_gap_i = not_gap.unsqueeze(1)
    not_gap_j = not_gap.unsqueeze(0)
    
    both_not_gap = not_gap_i & not_gap_j
    matches = (a3m_i == a3m_j) & both_not_gap
    
    # Count matches and totals
    num_matches = matches.sum(dim=2, dtype=torch.float32)
    num_positions = both_not_gap.sum(dim=2, dtype=torch.float32)
    
    # Compute identities
    identities = num_matches / num_positions.clamp(min=1)
    identities = torch.where(num_positions > 0, identities, torch.zeros_like(identities))
    
    # Diversity: average pairwise identity
    triu_mask = torch.triu(torch.ones(num_sequences, num_sequences, device=device), diagonal=1).bool()
    pairwise_identities = identities[triu_mask]
    diversity = float(1 - pairwise_identities.mean().cpu().item())
    
    # Neff: number of clusters at 80% identity
    num_neighbors = (identities >= 0.8).sum(dim=1).cpu().numpy()
    neff = float((1.0 / num_neighbors).sum())
    
    # Coverage: fraction of non-gap positions
    coverage = float((a3m != '-').mean())
    
    # Clean up GPU memory
    del a3m_tensor, not_gap, a3m_i, a3m_j, not_gap_i, not_gap_j, both_not_gap, matches
    torch.cuda.empty_cache()
    
    return {
        'M': num_sequences,
        'L': seq_len,
        'N_eff': neff,
        'diversity': diversity,
        'mean_coverage': coverage,
        'avg_length': avg_length
    }


def test_npz2msa():
    npz_path = '/mnt/ssd/public/sseqs/eval/runs_and_poses/boltz_results_t30_sseqs/processed/msa/7yha__1__1.A__1.E_1.F_1.G_1.H_0.npz'
    a3m_path = '/mnt/ssd/public/sseqs/eval/runs_and_poses/a3m/7yha__1__1.A__1.E_1.F_1.G_1.H.a3m'

    print(f"Testing NPZ to MSA conversion...")
    print(f"  NPZ: {npz_path}")
    print(f"  A3M: {a3m_path}")
    
    # Extract sequences from NPZ
    npz_seqs = npz2msa(npz_path)
    print(f"\n  NPZ file has {len(npz_seqs)} sequences")
    
    # Parse A3M file and remove lowercase (insertions)
    a3m_seqs = []
    with open(a3m_path, 'r') as f:
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    # Remove lowercase letters (insertions)
                    seq = ''.join(c for c in ''.join(current_seq) if not c.islower())
                    a3m_seqs.append(seq)
                    current_seq = []
            else:
                current_seq.append(line)
        if current_seq:
            seq = ''.join(c for c in ''.join(current_seq) if not c.islower())
            a3m_seqs.append(seq)
    
    print(f"  A3M file has {len(a3m_seqs)} sequences (lowercase removed)")
    
    # Compare sequences
    if len(npz_seqs) != len(a3m_seqs):
        print(f"\n✗ MISMATCH: Different number of sequences!")
        return False
    
    matches = 0
    mismatches = []
    for i in range(len(npz_seqs)):
        if npz_seqs[i] == a3m_seqs[i]:
            matches += 1
        else:
            if len(mismatches) < 3:  # Only store first 3 mismatches
                mismatches.append((i, npz_seqs[i][:50], a3m_seqs[i][:50]))
    
    if matches == len(npz_seqs):
        print(f"\n✓ ALL {matches} SEQUENCES MATCH!")
        return True
    else:
        print(f"\n✗ MISMATCH: {matches}/{len(npz_seqs)} sequences match")
        for i, npz, a3m in mismatches:
            print(f"  Sequence {i}:")
            print(f"    NPZ: {npz}...")
            print(f"    A3M: {a3m}...")
        return False

def test_neff():
    from icecream import ic 
    import os 
    path = '/mnt/ssd/public/sseqs/eval/runs_and_poses/low_homology_a3m_all/'
    files = sorted([a for a in os.listdir(path) if a[-4:] == '.a3m'])[:10]  # Just test first 10
    for fname in files: 
        pdb = fname[:4]
        p = f"{path}{fname}"
        max_lines = 2048

        # NEFFY-style with threshold 0.67 (excluding gaps from identity)
        us = neff(p, threshold=0.67, max_lines=max_lines)
        ic(pdb, round(us, 4))


# ============================================================================
# Fast vectorized PSSM computation  
# ============================================================================
def pssm_from_file(path, pca=1.0, pcb=1.5, bit_factor=8.0):
    """Fast PSSM: a3m/csv -> (query_str, PSSM[L,20] int8)"""
    AA = 'ARNDCQEGHILKMFPSTWYV'
    BG = np.array([.0787,.0512,.0448,.0537,.0151,.0404,.0633,.0747,.0226,
                   .0575,.0941,.0587,.0238,.0399,.0473,.0657,.0534,.0117,.0305,.0675], dtype=np.float32)
    
    # Parse file - keep as bytes, filter lowercase with numpy
    text = open(path, 'rb').read()
    if path.endswith('.csv'):
        lines = text.strip().split(b'\n')[1:]
        seqs_bytes = [l.split(b',')[1] for l in lines if b',' in l]
    else:
        lines = text.strip().split(b'\n')
        seqs_bytes = [lines[i+1] for i in range(0, len(lines)-1, 2) if lines[i].startswith(b'>')]
    if len(seqs_bytes) < 2: return None, None
    
    # Remove lowercase (97-122) using numpy - much faster than Python string ops
    query_raw = seqs_bytes[0]
    filtered = []
    for s in seqs_bytes:
        arr = np.frombuffer(s, dtype=np.uint8)
        filtered.append(bytes(arr[(arr < 97) | (arr > 122)]))
    seqs_bytes = filtered
    
    # Build matrix directly from bytes
    N, L = len(seqs_bytes), len(seqs_bytes[0])
    lut = np.full(256, -1, dtype=np.int8)
    for i, a in enumerate(AA): lut[ord(a)] = i
    mat = lut[np.frombuffer(b''.join(seqs_bytes), dtype=np.uint8)].reshape(N, L)
    valid = (mat >= 0)
    mat_safe = mat * valid  # faster than np.where
    
    # Counts using bincount (faster than add.at)
    flat_idx = (np.arange(L)[None, :] * 20 + mat_safe).ravel()  # pos*20 + aa
    flat_valid = valid.ravel().astype(np.float32)
    counts = np.bincount(flat_idx, weights=flat_valid, minlength=L*20).reshape(L, 20).astype(np.float32)
    distinct = np.maximum((counts > 0).sum(1), 1.0)
    
    # Henikoff weights - avoid where, use multiply
    nres = (valid.sum(1) + 30.0).astype(np.float32)[:, None]
    pos_idx = np.broadcast_to(np.arange(L), (N, L))
    seq_aa_counts = counts[pos_idx.ravel(), mat_safe.ravel()].reshape(N, L)
    seq_aa_counts = seq_aa_counts + (1.0 - valid)  # set invalid to 1 to avoid div by 0
    w = (valid / (seq_aa_counts * distinct * nres)).sum(1)
    w = ((w + 1e-6) / (w.sum() + 1e-6)).astype(np.float32)
    
    # Weighted frequencies using bincount
    flat_idx2 = (np.arange(L)[None, :] * 20 + mat_safe).ravel()
    flat_w = (w[:, None] * valid).ravel()
    freq = np.bincount(flat_idx2, weights=flat_w, minlength=L*20).reshape(L, 20).astype(np.float32)
    freq /= np.maximum(freq.sum(1, keepdims=True), 1e-10)
    
    # Neff + pseudocounts + log-odds
    ent = -np.sum(np.where(freq > 1e-10, freq * np.log2(np.maximum(freq, 1e-10)), 0), 1)
    tau = np.minimum(1.0, pca / (1.0 + (2.0 ** ent) / pcb))[:, None]
    profile = (1 - tau) * freq + tau * BG
    pssm = np.clip(np.round(bit_factor * np.log2(np.maximum(profile, 1e-10) / BG)), -128, 127).astype(np.int8)
    
    return filtered[0].decode('latin1'), pssm


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='MSA utilities')
    parser.add_argument('-test_npz2msa', action='store_true') 
    parser.add_argument('-test_neff', action='store_true') 
    parser.add_argument('-test_pssm', action='store_true', help='Test PSSM on a3m/csv file')
    args = parser.parse_args()
    
    if args.test_npz2msa: test_npz2msa()
    if args.test_neff: test_neff()
    if args.test_pssm:
        import time
        t0 = time.time()
        q, pssm = pssm_from_file()
        t = time.time() - t0
        if pssm is not None:
            print(f"Query: {len(q)} aa, PSSM: {pssm.shape}, range: [{pssm.min()}, {pssm.max()}], time: {t*1000:.1f}ms")
