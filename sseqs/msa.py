import os 
import sys 
import hashlib
import pickle 
import numpy as np 
import math 
from torch.utils.cpp_extension import load
import torch as th
from tqdm import tqdm
import psutil

# CHUNKS = number of GB to load into RAM. 
chunks = int(os.environ.get('CHUNKS', 32))
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

    # Calculate chunks to load (min of requested and available in DB)
    total_db_size = sum(os.path.getsize(f) for f in xbit_files)
    total_chunks_in_db = math.ceil(total_db_size / 1_000_000_000)
    num_chunks = min(chunks, total_chunks_in_db)

    # Check available RAM
    avail_gb_ram = psutil.virtual_memory().available / 1_000_000_000
    needed_gb_ram = num_chunks + 10
    if avail_gb_ram < needed_gb_ram: 
        print(f"Need {needed_gb_ram}GB RAM for {num_chunks}GB database, only {avail_gb_ram:.1f}GB available.")
        exit()

    # Prepare pinned RAM to enable fast RAM->VRAM transfer.
    pinned_ram = [th.empty((1_000_000_000), dtype=th.uint8, pin_memory=True) for i in tqdm(range(num_chunks), desc=f"Preparing {num_chunks}GB of pinned RAM")]

    # Read all chunks from files into pinned RAM
    desc = f"Reading {num_chunks}GB MSA data from {len(xbit_files)} file(s) into pinned RAM"
    pbar = tqdm(total=num_chunks, desc=desc, unit='GB')
    chunk_idx = 0
    for xbit_file in sorted(xbit_files):
        with open(xbit_file, "rb") as fh:
            file_size = os.path.getsize(xbit_file)
            file_chunks = file_size // (1000 * 1000 * 1000)
            for file_chunk in range(file_chunks):
                if chunk_idx >= num_chunks:
                    break
                tensor = pinned_ram[chunk_idx]
                fh.seek(file_chunk * 1000 * 1000 * 1000)
                fh.readinto(tensor.numpy())
                pbar.update(1)
                chunk_idx += 1
    pbar.close()

    sys.argv.append('-loaded')

from sseqs.to_a3m import to_a3m

import threading
_last={}

def msa(query, filenames=None, num_msas: int=8192, save_raw: str='', 
        verbose: bool=True, 
        bs: int=100_000_000, 
        sw_threshold: int=70,  # keep hits with SW score > threshold
        top_k: int=100_000,  # keep top K hits
        gpu: int = 0, 
        blocking: bool = True): 
    global pinned_ram

    # Support multiple GPUs reading MSA from same RAM asynchronously (non-blocking mode)
    if not blocking: 
        print(gpu)
        t=_last.get(gpu)
        if t is not None and t.is_alive(): t.join()
        def _target():
            msa(
                query=query,
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

    with th.cuda.device(gpu):
        queries = [query]
        if type(filenames) == str: filenames = [filenames]

        if filenames == None: 
            hashes = [hashlib.sha256(query.encode()).hexdigest()[:8] for query in queries]
            filenames = [hash + '.a3m' for hash in hashes]

        aas = "ARNDCQEGHILKMFPSTWYV"                      # 20 standard AAs
        AA_TO_POS = {aa: i for i, aa in enumerate(aas)}
        AA_TO_POS.update({"X": 20, "-": 20, "*": 20})      # gap / unknown

        max_sw, matches = 0, 0
        dbhits = []
        scores = []
        AA_MAP = {aa: idx for idx, aa in enumerate("ARNDCQEGHILKMFPSTWYV")}

        query_enc = th.tensor([AA_MAP[aa] for aa in query if aa in AA_MAP], dtype=th.uint8, device=f'cuda:{gpu}')

        arange = th.arange(int(bs*1.6), device=f'cuda:{gpu}', dtype=th.int32)

        # Create separate streams to overlap RAM→VRAM copying and compute
        copy_stream = th.cuda.Stream(device=f'cuda:{gpu}')
        compute_stream = th.cuda.Stream(device=f'cuda:{gpu}')
        
        current_chunk = th.zeros_like(pinned_ram[0], device=f'cuda:{gpu}')
        next_chunk = th.empty_like(pinned_ram[0], device=f'cuda:{gpu}')
        
        # Track cumulative statistics
        total_seqs_processed = 0
        total_gcups_processed = 0
        total_hits = 0
        last_mode_len = 0

        total_chunks = len(pinned_ram)
        
        if verbose: pbar = tqdm(range(total_chunks), unit='GB')
        else: pbar = range(total_chunks)

        for i in pbar:
            a_pin = pinned_ram[i]

            # Double-buffer: copy next chunk while processing current
            if i != 0:
                copy_stream.synchronize()
                current_chunk, next_chunk = next_chunk, current_chunk
            else:
                with th.cuda.stream(copy_stream):
                    current_chunk.copy_(a_pin, non_blocking=True)
                copy_stream.synchronize()
            
            # Start async copy of next chunk
            if i + 1 < total_chunks:
                with th.cuda.stream(copy_stream):
                    next_chunk.copy_(pinned_ram[i + 1], non_blocking=True)

            # Compute
            with th.cuda.stream(compute_stream):
                results = [] 
        
                assert current_chunk.numel() % bs == 0, "current_chunk.numel() must be divisible by bs"

                for j in range(0, math.ceil(current_chunk.numel()/bs), 1):
                    ascii = xbit.decompress(current_chunk, int(bs*8/5), j*bs)
                    delim_pos = arange[ascii==64].flatten()

                    lengths = delim_pos[1:] - delim_pos[:-1] - 1

                    valid_lengths = lengths[lengths > 0]
                    if valid_lengths.numel() == 0: continue
                    
                    max_len = th.max(valid_lengths)
                    
                    good_lengths_mask = (lengths > 0)
                    d_indices = th.nonzero(good_lengths_mask, as_tuple=True)[0]
                    
                    if d_indices.numel() == 0: continue
                    
                    # Skip first and last sequence to ensure clean boundaries
                    if d_indices.numel() <= 2: continue
                    d_indices = d_indices[1:-1]
                    
                    # Round down to multiple of SEQS_PER_BLOCK
                    SEQS_PER_BLOCK = 64
                    num_keep = (d_indices.numel() // SEQS_PER_BLOCK) * SEQS_PER_BLOCK
                    if num_keep == 0: continue
                    d_indices = d_indices[:num_keep]
                    
                    d = delim_pos[d_indices + 1]
                    first_idx = d_indices[0] if d_indices.numel() > 0 else 0
                    first_delim = delim_pos[first_idx].view(1)
                    
                    starts = th.cat([first_delim.view(1), d]).to(th.int32).contiguous()
                    good_idx = th.arange(1, d.numel() + 1, dtype=th.int32, device=f'cuda:{gpu}')
                    
                    sw_scores = sw.sw_cuda_affine(query_enc, good_idx, ascii, starts, max_len, 11, 1)
                    
                    batch_gcups = float(good_idx.numel()) * float(max_len.item()) * len(query) / 1e9
                    total_seqs_processed += good_idx.numel()
                    total_gcups_processed += batch_gcups
                    last_mode_len = max_len
                                
                    mask = sw_scores > sw_threshold
                    if not mask.any(): continue
                    
                    sw_scores = sw_scores[mask]
                    topk_idx = th.nonzero(mask, as_tuple=True)[0]
                    
                    total_hits += topk_idx.numel()
                    
                    topk_good_idx = good_idx[topk_idx]
                    seq_starts = starts[topk_good_idx - 1] + 1
                    seq_ends = starts[topk_good_idx]
                    seq_lengths = seq_ends - seq_starts
                    
                    seq_starts_i64 = seq_starts.to(th.int64)
                    idx = seq_starts_i64[:, None] + th.arange(max_len, device=f'cuda:{gpu}')[None, :]
                    seqs = ascii.expand(idx.size(0), -1).gather(1, idx)
                    
                    results.append((j, i, sw_scores, seq_lengths, seqs, max_len))
            
            compute_stream.synchronize()

            for j, i, sw_scores, seq_lengths, current_batch_seqs, max_len_val in results: 
                B = current_batch_seqs.shape[1]
                buf_u8 = current_batch_seqs.cpu().numpy()
                n, B    = buf_u8.shape
                rows_bytes = buf_u8.view(f'S{B}').ravel()
                rows_unicode = rows_bytes.astype(f'U{B}')
                lengths = seq_lengths.cpu().numpy().tolist()
                current_batch_seqs = [s[:l] for s, l in zip(rows_unicode, lengths)]

                max_sw_item = th.max(sw_scores).item()
                max_sw = max(max_sw, max_sw_item)
                matches += len(current_batch_seqs)

                dbhits.append(current_batch_seqs)
                scores.append(sw_scores.cpu().numpy())

            # Prune to top_k if we've accumulated too many hits
            if top_k and total_hits > top_k * 2:
                scores_cat = np.concatenate(scores)
                dbhits_cat = np.concatenate(dbhits)
                idx = np.argpartition(scores_cat, -top_k)[-top_k:]
                scores = [scores_cat[idx]]
                dbhits = [dbhits_cat[idx]]
                total_hits = top_k

            if verbose:
                dblen = int(last_mode_len.item()) if isinstance(last_mode_len, th.Tensor) else last_mode_len
                pbar.set_description(
                    f"qlen={len(query)} dblen={dblen} seqs={total_seqs_processed/1e6:.1f}M hits={total_hits/1e3:.1f}K GCUPs={total_gcups_processed:.0f}"
                )
        
        query, dbhits, scores = queries[0], np.concatenate(dbhits), np.concatenate(scores)

        # Final top-k selection
        if top_k and len(scores) > top_k:
            idx = np.argpartition(scores, -top_k)[-top_k:]
            idx = idx[np.argsort(scores[idx])[::-1]]
            scores = scores[idx]
            dbhits = dbhits[idx]

        if save_raw != '': 
            if verbose: print(f"Saving to `{save_raw}`.")
            pickle.dump([query, [hit.encode('ascii') for hit in dbhits], scores], open(save_raw, 'wb'))

        idx = np.argsort(scores)[::-1] 
        _seq_matches = dbhits[idx][:num_msas]

        if filenames[0] is not None: 
            if verbose: 
                print(f"Saving result to {filenames[0]}")
            to_a3m(query, _seq_matches, filename=filenames[0]) 


if __name__ == "__main__":
    msa('ADAM'*30, "adam.a3m")
    msa('ADA'*30, "ada.a3m")
