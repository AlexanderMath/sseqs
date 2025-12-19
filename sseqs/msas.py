# like `msa.py` but for running batch of queries => no loading to ram just load to disk. 
import os 
import sys 
import pickle
from typing import List 
import numpy as np 
import math 
from torch.utils.cpp_extension import load
import torch as th
from tqdm import tqdm

# Loads `CHUNKS` GB of MSA database to process (streams from disk)
chunks = int(os.environ.get('CHUNKS', 4))
assert chunks > 0

# Look for database. 
xbit_path = os.environ.get('DBPATH', 'uniref_bfd_mgy_cf.xbit')
package_dir = os.path.dirname(os.path.abspath(__file__))

# @alex:  supports cache loading in ipython notebook
#         reach out if you can find nicer way to achieve this. 
if sys.argv[-1] != '-loaded': 

    # Auto-detect GPU and compile only for that architecture
    if th.cuda.is_available():
        major, minor = th.cuda.get_device_capability()
        os.environ["TORCH_CUDA_ARCH_LIST"] = f"{major}.{minor}"

    # Get list of xbit files: if path is a directory, load all .xbit files; otherwise just the file
    if os.path.isdir(xbit_path):
        xbit_files = sorted([os.path.join(xbit_path, f) for f in os.listdir(xbit_path) if f.endswith('.xbit')])
        if not xbit_files:
            print(f"No .xbit files found in directory `{xbit_path}`")
            exit()
    else:
        if not os.path.exists(xbit_path):
            print(f"Didn't find MSA database at `{xbit_path}`")
            print("\twget https://foldify.org/uniref_bfd_mgy_cf.xbit")
            print("\texport DBPATH=$PWD/uniref_bfd_mgy_cf.xbit")
            exit()
        xbit_files = [xbit_path]

    # Load CUDA kernels. 
    verbose = True
    xbit = load("xbit", sources=[f"{package_dir}/xbit_bind.cpp",   f"{package_dir}/xbit.cu"],   extra_cuda_cflags=["-O3", "--use_fast_math"], extra_cflags=["-O3"], verbose=verbose)
    sw   = load("sw",   sources=[f"{package_dir}/sw_bind.cpp",     f"{package_dir}/sw.cu"],     extra_cuda_cflags=["-O3", "--use_fast_math"], extra_cflags=["-O3"], verbose=verbose)

    # Batch mode: no RAM buffer, stream directly from disk
    print(f"Batch mode: will stream {chunks}GB from disk ({len(xbit_files)} file(s))")

    sys.argv.append('-loaded')

from sseqs.to_a3m import to_a3m

import threading
_last={}

def msas(queries: List[str], filenames: List[str], num_msas: int=8192, save_raw: List[str] = None, 
        verbose: bool=True, 
        bs: int=100_000_000, 
        sw_threshold: int=70,  # keep hits with SW score > threshold (or None for top_k mode)
        top_k: int=100_000,  # keep top K hits (streaming mode)
        gpu: int = 0, 
        blocking: bool = True): 

    # Support 4 gpu reading MSA from same RAM asynchronously = non-blocking
    if not blocking: 
        print(gpu)
        t=_last.get(gpu)
        if t is not None and t.is_alive(): t.join()
        # Capture all variables explicitly for the closure
        def _target():
            msas(
                queries=queries,
                filenames=filenames,
                num_msas=num_msas,
                save_raw=save_raw,
                verbose=verbose,
                bs=bs,
                sw_threshold=sw_threshold,
                top_k=top_k,
                gpu=gpu,
                blocking=True)
        nt=threading.Thread(target=_target)
        _last[gpu]=nt
        nt.start()
        return nt

    # Support different gpu in blocking/sequential mode. 
    with th.cuda.device(gpu):
        # Validate inputs
        assert isinstance(queries, list) and len(queries) > 0, "queries must be a non-empty list"
        assert isinstance(filenames, list) and len(filenames) == len(queries), f"filenames length must match queries: {len(filenames)} vs {len(queries)}"
        if save_raw is not None: assert isinstance(save_raw, list) and len(save_raw) == len(queries), f"save_raw length must match queries: {len(save_raw)} vs {len(queries)}"

        aas = "ARNDCQEGHILKMFPSTWYV"                      # 20 standard AAs
        AA_TO_POS = {aa: i for i, aa in enumerate(aas)}
        AA_TO_POS.update({"X": 20, "-": 20, "*": 20})      # gap / unknown

        # Initialize per-query storage
        all_query_data = []
        for query_idx, query in enumerate(queries):
            AA_MAP = {aa: idx for idx, aa in enumerate("ARNDCQEGHILKMFPSTWYV")}
            query_enc = th.tensor([AA_MAP[aa] for aa in query if aa in AA_MAP], dtype=th.uint8, device=f'cuda:{gpu}')
            all_query_data.append({
                'query': query,
                'query_enc': query_enc,
                'dbhits': [],
                'scores': [],
                'max_sw': 0,
                'matches': 0,
                'total_seqs_processed': 0,
                'total_gcups_processed': 0,
                'total_hits': 0,
                'last_mode_len': 0
            })

        arange = th.arange(int(bs*1.6), device=f'cuda:{gpu}', dtype=th.int32)

        # Create separate streams to overlap disk I/O and compute
        copy_stream = th.cuda.Stream(device=f'cuda:{gpu}')
        compute_stream = th.cuda.Stream(device=f'cuda:{gpu}')
        
        # Calculate total chunks to process from disk
        total_db_size = sum(os.path.getsize(f) for f in xbit_files)
        total_chunks_in_db = math.ceil(total_db_size / 1_000_000_000)
        total_chunks_to_process = min(chunks, total_chunks_in_db)
        
        # Allocate buffers: one pinned for disk reads, two on GPU for double buffering
        temp_pinned = th.empty((1_000_000_000), dtype=th.uint8, pin_memory=True)
        current_chunk = th.zeros((1_000_000_000), dtype=th.uint8, device=f'cuda:{gpu}')
        next_chunk = th.empty((1_000_000_000), dtype=th.uint8, device=f'cuda:{gpu}')
        
        # Helper to load chunk from disk into pinned buffer
        def load_chunk_from_disk(chunk_id, buffer):
            target_byte = chunk_id * 1_000_000_000
            cumulative = 0
            for fname in sorted(xbit_files):
                fsize = os.path.getsize(fname)
                if target_byte < cumulative + fsize:
                    file_offset = target_byte - cumulative
                    with open(fname, 'rb') as f:
                        f.seek(file_offset)
                        f.readinto(buffer.numpy())
                    return
                cumulative += fsize

        # Remove output files
        for filename in filenames:
            if os.path.exists(filename): os.remove(filename)

        if verbose: pbar = tqdm(range(total_chunks_to_process), unit='GB')
        else: pbar = range(total_chunks_to_process)

        for i in pbar:
            # Load chunk from disk into pinned buffer
            load_chunk_from_disk(i, temp_pinned)
            
            # Copy to GPU (with double buffering)
            if i != 0:
                # Wait for previous async copy to finish, then swap
                copy_stream.synchronize()
                current_chunk, next_chunk = next_chunk, current_chunk
            else:
                # First iteration: copy data into pre-allocated buffer
                with th.cuda.stream(copy_stream):
                    current_chunk.copy_(temp_pinned, non_blocking=True)
                copy_stream.synchronize()
            
            # Start NEXT async copy immediately (load from disk while computing)
            if i + 1 < total_chunks_to_process:
                # Load next chunk from disk in background (this could be threaded for more overlap)
                load_chunk_from_disk(i + 1, temp_pinned)
                with th.cuda.stream(copy_stream):
                    next_chunk.copy_(temp_pinned, non_blocking=True)

            # Do the compute bit - FOR ALL QUERIES
            with th.cuda.stream(compute_stream):

                # Process all queries against this chunk
                for query_idx, qdata in enumerate(all_query_data):
                    print('\r', query_idx,  end='', flush=True)
                    results = [] 
                    query_enc = qdata['query_enc']
                    query = qdata['query']
            
                    assert current_chunk.numel() % bs == 0, "current_chunk.numel() must be divisible by bs"

                    for j in range(0, math.ceil(current_chunk.numel()/bs), 1):
                        ascii = xbit.decompress(current_chunk, int(bs*8/5), j*bs)
                        delim_pos = arange[ascii==64].flatten()

                        # Keep only sequences with mode length
                        lengths = delim_pos[1:] - delim_pos[:-1] - 1

                        # Filter to valid lengths
                        valid_lengths = lengths[lengths > 0]
                        if valid_lengths.numel() == 0: continue
                        
                        # Sequences are sorted by length, so just use max for dispatch
                        max_len = th.max(valid_lengths)
                        
                        # Keep ALL valid sequences (they're already close in length since sorted)
                        good_lengths_mask = (lengths > 0)
                        d_indices = th.nonzero(good_lengths_mask, as_tuple=True)[0]
                        
                        if d_indices.numel() == 0: continue
                        
                        # Skip first and last sequence to ensure clean boundaries
                        if d_indices.numel() <= 2: continue
                        d_indices = d_indices[1:-1]
                        
                        # Round down to multiple of SEQS_PER_BLOCK to avoid out-of-bounds writes
                        SEQS_PER_BLOCK = 64  # Conservative value (kernel uses 16-64 depending on length)
                        num_keep = (d_indices.numel() // SEQS_PER_BLOCK) * SEQS_PER_BLOCK
                        if num_keep == 0: continue
                        d_indices = d_indices[:num_keep]
                        
                        # d_indices are indices into lengths array, where lengths[i] = sequence between delim_pos[i] and delim_pos[i+1]
                        # So delim_pos[d_indices[i]] is delimiter BEFORE sequence, delim_pos[d_indices[i]+1] is delimiter AFTER
                        # We need delimiters AFTER sequences for starts array
                        d = delim_pos[d_indices + 1]  # Delimiters AFTER each kept sequence
                        # Get delimiter BEFORE first sequence (d_indices[0] is the index into lengths array)
                        # Use clamp to avoid sync: if d_indices[0] < 0, clamp it to 0 and we'll get wrong value but mask it later
                        first_idx = d_indices[0] if d_indices.numel() > 0 else 0
                        first_delim = delim_pos[first_idx].view(1)  # No if statement - stays on GPU!
                        
                        # starts array: [first_delim, d[0], d[1], ..., d[N-1]]
                        # d[-1] is already the delimiter after the last sequence, so we don't need a separate last_delim
                        starts = th.cat([first_delim.view(1), d]).to(th.int32).contiguous()
                        good_idx = th.arange(1, d.numel() + 1, dtype=th.int32, device=f'cuda:{gpu}')
                        
                        # Kernel call (pass max_len for dispatch - kernel reads actual lengths from starts array)
                        sw_scores = sw.sw_cuda_affine(query_enc, good_idx, ascii, starts, max_len, 11, 1)
                        
                        # Update cumulative stats (convert to Python types to avoid overflow)
                        batch_gcups = float(good_idx.numel()) * float(max_len.item()) * len(query) / 1e9
                        qdata['total_seqs_processed'] += good_idx.numel()
                        qdata['total_gcups_processed'] += batch_gcups
                        qdata['last_mode_len'] = max_len  # Store for progress bar (still a tensor, no sync)
                                    
                        # Filter by SW threshold
                        mask = sw_scores > sw_threshold
                        if not mask.any(): continue
                        
                        sw_scores = sw_scores[mask]
                        topk_idx = th.nonzero(mask, as_tuple=True)[0]
                        
                        # Track hits
                        qdata['total_hits'] += topk_idx.numel()
                        
                        # Extract sequences - now with VARIABLE lengths!
                        # Get the actual start/end positions for each sequence from good_idx
                        topk_good_idx = good_idx[topk_idx]  # good_idx indices for top-k sequences
                        seq_starts = starts[topk_good_idx - 1] + 1  # Start of each sequence
                        seq_ends = starts[topk_good_idx]  # End of each sequence (delimiter position)
                        seq_lengths = seq_ends - seq_starts  # Actual length of each sequence
                        
                        # Extract using max_len as the buffer size (some will be padded)
                        seq_starts_i64 = seq_starts.to(th.int64)
                        idx = seq_starts_i64[:, None] + th.arange(max_len, device=f'cuda:{gpu}')[None, :]
                        seqs = ascii.expand(idx.size(0), -1).gather(1, idx)
                        
                        # Store sequences with their actual lengths (not just max_len)
                        results.append((j, i, sw_scores, seq_lengths, seqs, max_len))
                    
                    # Process results for this query
                    for j, i, sw_scores, seq_lengths, current_batch_seqs, max_len_val in results: 
                        B = current_batch_seqs.shape[1]
                        buf_u8 = current_batch_seqs.cpu().numpy()      # (n, B)  dtype=uint8
                        n, B    = buf_u8.shape
                        rows_bytes = buf_u8.view(f'S{B}').ravel()      # (n,)   dtype='|S<B>'
                        rows_unicode = rows_bytes.astype(f'U{B}')      # (n,)   dtype='<U<B>'
                        # Now we can safely sync to get the actual lengths for each sequence
                        lengths = seq_lengths.cpu().numpy().tolist()  # Variable lengths!
                        current_batch_seqs = [s[:l] for s, l in zip(rows_unicode, lengths)]

                        # --- OUT ---------------------------------------------
                        max_sw_item = th.max(sw_scores).item()  # Single sync point per batch
                        qdata['max_sw'] = max(qdata['max_sw'], max_sw_item)
                        qdata['matches'] += len(current_batch_seqs)

                        qdata['dbhits'].append(current_batch_seqs)
                        qdata['scores'].append(sw_scores.cpu().numpy())

                    # Prune to top_k if we've accumulated too many hits for this query
                    if top_k and qdata['total_hits'] > top_k * 2:
                        scores_cat = np.concatenate(qdata['scores'])
                        dbhits_cat = np.concatenate(qdata['dbhits'])
                        idx = np.argpartition(scores_cat, -top_k)[-top_k:]
                        qdata['scores'] = [scores_cat[idx]]
                        qdata['dbhits'] = [dbhits_cat[idx]]
                        qdata['total_hits'] = top_k

            
            # Wait for compute to finish before processing results
            compute_stream.synchronize()

            # Update progress bar with stats from first query
            if verbose and len(all_query_data) > 0:
                qdata = all_query_data[0]
                dblen = int(qdata['last_mode_len'].item()) if isinstance(qdata['last_mode_len'], th.Tensor) else qdata['last_mode_len']
                pbar.set_description(
                    f"{len(queries)}Q qlen={len(qdata['query'])} dblen={dblen} seqs={qdata['total_seqs_processed']/1e6:.1f}M hits={qdata['total_hits']/1e3:.1f}K GCUPs={qdata['total_gcups_processed']:.0f}"
                )

        # Finalize all queries
        for query_idx, qdata in enumerate(all_query_data):
            query = qdata['query']
            dbhits = np.concatenate(qdata['dbhits']) if qdata['dbhits'] else np.array([])
            scores = np.concatenate(qdata['scores']) if qdata['scores'] else np.array([])

            # Final top-k selection
            if top_k and len(scores) > top_k:
                idx = np.argpartition(scores, -top_k)[-top_k:]
                idx = idx[np.argsort(scores[idx])[::-1]]  # Sort the top-k
                scores = scores[idx]
                dbhits = dbhits[idx]
            if save_raw is not None and save_raw[query_idx]: 
                if verbose: print(f"Saving to `{save_raw[query_idx]}`.")
                pickle.dump([query, [hit.encode('ascii') for hit in dbhits], scores], open(save_raw[query_idx], 'wb'))

            idx = np.argsort(scores)[::-1] 
            _seq_matches = dbhits[idx][:num_msas]

            if filenames[query_idx] is not None: 
                if verbose: 
                    print(f"Saving result to {filenames[query_idx]}")
                to_a3m(query, _seq_matches, filename=filenames[query_idx]) 



def find_sequences(targets: List[str], gpu: int = 0, verbose: bool = True) -> dict:
    """Check if exact sequences exist in the database. Returns dict {seq: True/False}."""
    
    # Pre-encode targets as bytes for fast set lookup
    targets_bytes = {seq.encode('ascii'): seq for seq in targets}
    found = {seq: False for seq in targets}
    remaining = len(targets)
    
    total_db_size = sum(os.path.getsize(f) for f in xbit_files)
    total_chunks = math.ceil(total_db_size / 1_000_000_000)
    
    temp_pinned = th.empty((1_000_000_000), dtype=th.uint8, pin_memory=True)
    current_chunk = th.zeros((1_000_000_000), dtype=th.uint8, device=f'cuda:{gpu}')
    
    def load_chunk_from_disk(chunk_id, buffer):
        target_byte = chunk_id * 1_000_000_000
        cumulative = 0
        for fname in sorted(xbit_files):
            fsize = os.path.getsize(fname)
            if target_byte < cumulative + fsize:
                file_offset = target_byte - cumulative
                with open(fname, 'rb') as f:
                    f.seek(file_offset)
                    f.readinto(buffer.numpy())
                return
            cumulative += fsize
    
    pbar = tqdm(range(total_chunks), desc="Searching") if verbose else range(total_chunks)
    
    for i in pbar:
        load_chunk_from_disk(i, temp_pinned)
        current_chunk.copy_(temp_pinned)
        
        # Decompress full chunk at once
        ascii_data = xbit.decompress(current_chunk, int(1_000_000_000 * 8 / 5), 0)
        all_bytes = ascii_data.cpu().numpy().tobytes()
        
        # Split by '@' once (C-level), build set of bytes
        db_seqs = set(all_bytes.split(b'@'))
        
        # O(1) lookup per target
        for seq_bytes, seq in targets_bytes.items():
            if not found[seq] and seq_bytes in db_seqs:
                found[seq] = True
                remaining -= 1
                if verbose:
                    print(f"\nFound ({len(targets)-remaining}/{len(targets)}): {seq[:50]}...")
        
        if remaining == 0:
            return found
        
        if verbose:
            pbar.set_postfix(found=f"{len(targets)-remaining}/{len(targets)}")
    
    return found


if __name__ == "__main__":
    msas(['ADAM'*30], ["adam.a3m"])
    msas(['ADA'*30], ["ada.a3m"])