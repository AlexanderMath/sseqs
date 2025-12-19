from sseqs.sw import sw_affine_backtrack
from tqdm import tqdm 
import numpy as np 

def to_a3m(query: str, seqs: list[str], filename: str, gap_open: int = 11, gap_extend: int = 1, batch_size: int = 1024):
    # todo: try with gap_open=10 ?
    """
    Write an A3M-formatted alignment where each target row is padded with
    leading/trailing gaps so that it covers the full query length.
    """
    aligned = []
    pbar = tqdm(range(0, len(seqs), batch_size))
    for start in pbar:
        stop = start + batch_size
        aligned = aligned + sw_affine_backtrack(query, seqs[start: stop], gap_open=gap_open, gap_extend=gap_extend)[-1]
        pbar.set_description(f'dblen={len(seqs[start])}')
    # merge into one list, then sort by score
    scores = [a['score'] for a in aligned]
    aligned = [aligned[i] for i in np.argsort(scores)[::-1]]


    # -------------------------------------------------------------------------
    # helper that converts ONE local alignment into a full-length A3M row
    # -------------------------------------------------------------------------
    def _build_target_row(align_res: dict) -> str:
        q_loc          = align_res["q_aligned"]     # with gaps
        t_loc          = align_res["t_aligned"]     # with gaps

        # ------- (2) gap padding BEFORE the alignment ------------------------
        q_compact = q_loc.replace('-', '')
        q_start_in_query = query.find(q_compact)            # 0-based
        if q_start_in_query == -1:
            raise ValueError("Cannot locate aligned fragment in query.")
        gap_prefix = '-' * q_start_in_query                 # columns 0 … q_start-1

        # ------- (3) the aligned block itself --------------------------------
        #        upper-case  = residues matching a query column
        #        lower-case  = insertions inside the block
        middle = []
        for q_char, t_char in zip(q_loc, t_loc):
            if t_char == '-':           # deletion in target
                middle.append('-')
            else:
                middle.append(t_char.upper() if q_char != '-' else t_char.lower())
        middle = ''.join(middle)

        # ------- (4) gap padding AFTER the alignment -------------------------
        gap_suffix = '-' * (len(query) - (q_start_in_query + len(q_compact)))

        # ------- concatenate everything --------------------------------------
        return gap_prefix + middle + gap_suffix
    # -------------------------------------------------------------------------

    #with open(filename, "a") as f:
    a3m = ""
    a3m += f">original_query\n{query}\n"
    seen = {}
    for i, ares in enumerate(aligned):
        if ares["score"] == 0 or not ares["q_aligned"]: continue # skip 

        row = _build_target_row(ares)
        if "@" in row: continue 

        # problem:  .a3m has multiple identical rows
        # cause:    different proteins have the same best alignment match 
        # solution: only add we didn't prev see 
        # @alex:    in some cases there'll be multiple identical loss alignments, so
        #           could give structure model more info given another one (or maybe 1% lower loss alignment)
        if row in seen: continue 
        seen[row] = 1
        a3m += f">target_{i} score={ares['score']}\n"
        a3m += row + "\n"

    open(filename, "w").write(a3m)