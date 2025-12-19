/*
    Linear Gapped Smith-Waterman using SIMD half2. 

        sw_affine(query, db[i], 11, 1) <= sw_linear(query, db[i], 1)

    Pre-filter provably doesn't throw away bad one. 
    Using uint8 provably does the same up to sw_score = 128 
    Use only positive blosum scores no negatives => fit registers and remove one max. 
    blosum only has 41 positive values in triagnular form 

*/
template<int BLOCK_DIM_X, int SUB_WARP_SIZE, int TPT>
__global__ void sw_kernel_linear( const unsigned char* query_seq_indices, int query_seq_len, int* good_idx, const uint8_t* ascii,
    int* starts, // call length earlier -- refactor to call starts, it contains the starting position of all sequences!
    int num_db_seqs, int* out_scores, int gap_penalty) {
    
#if COMPILE_LINEAR
    
    // DOUBLE-PACKING with SUB-WARPS: Each sub-warp processes 2 database sequences
    // half2.x = sequence 0, half2.y = sequence 1 (different sequences, same query position!)
    constexpr int WARP_SIZE = 32;
    
    const int warp_id = threadIdx.x / WARP_SIZE;
    const int lane = threadIdx.x % WARP_SIZE;
    const int sub_warp_id = lane / SUB_WARP_SIZE;  // Which sub-warp within the warp (0-3)
    const int sub_lane = lane % SUB_WARP_SIZE;      // Position within sub-warp (0-7)
    
    constexpr int WARPS_PER_BLOCK = BLOCK_DIM_X / WARP_SIZE;
    constexpr int SUB_WARPS_PER_WARP = WARP_SIZE / SUB_WARP_SIZE;  // 4 sub-warps per warp
    constexpr int SEQS_PER_SUB_WARP = 2;  // Double-packing: 2 sequences per sub-warp
    constexpr int SEQS_PER_WARP = SUB_WARPS_PER_WARP * SEQS_PER_SUB_WARP;  // 8 sequences per warp
    constexpr int SEQS_PER_BLOCK = WARPS_PER_BLOCK * SEQS_PER_WARP;
    
    // Each sub-warp processes 2 sequences
    const int seq0_idx = blockIdx.x * SEQS_PER_BLOCK + warp_id * SEQS_PER_WARP + sub_warp_id * 2;
    const int seq1_idx = blockIdx.x * SEQS_PER_BLOCK + warp_id * SEQS_PER_WARP + sub_warp_id * 2 + 1;
    
    // Load 3D combined BLOSUM into shared memory
    __shared__ half2 s_blosum_combined[NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA];
    
    // Cooperative load - ALL threads must participate before __syncthreads()
    for (int i = threadIdx.x; i < NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA; i += BLOCK_DIM_X) {
        s_blosum_combined[i] = blosum62_combined_cuda_global[i];
    }
    __syncthreads();
    
    // NOW we can do bounds checks
    bool seq_valid = (seq0_idx < num_db_seqs);
    
    // Check if seq1 is valid (for odd number of sequences)
    bool seq1_valid = (seq1_idx < num_db_seqs);
    
    // Get sequence 0 info (assume valid - host should pad if needed)
    int start0 = seq_valid ? starts[good_idx[seq0_idx]-1] + 1 : 0;
    int stop0 = seq_valid ? starts[good_idx[seq0_idx]] + 1 : 1;
    int length0 = stop0 - start0 - 1;
    const uint8_t* ascii0 = ascii + start0;
    
    // Get sequence 1 info (handle case where seq1 doesn't exist for odd counts)
    int start1, stop1, length1;
    const uint8_t* ascii1;
    if (seq1_valid) {
        start1 = starts[good_idx[seq1_idx]-1] + 1;
        stop1 = starts[good_idx[seq1_idx]] + 1;
        length1 = stop1 - start1 - 1;
        ascii1 = ascii + start1;
    } else {
        // Dummy values for seq1 (will be masked out)
        start1 = start0;
        stop1 = stop0;
        length1 = 0;  // Set to 0 so all positions are out-of-bounds
        ascii1 = ascii0;
    }
    
    // TRANSPOSED: Each thread handles TPT consecutive DB positions
    const int db_begin = sub_lane * TPT;  // Use sub_lane for sub-warp processing

    // Mark invalid if this thread's entire range is out of bounds (NO early return!)
    if (!seq_valid || (db_begin >= length0 && db_begin >= length1)) {
        length0 = 0;
        length1 = 0;
    }

    // Compute sub-warp mask for shuffle operations (only communicate within sub-warp)
    const unsigned int sub_warp_mask = ((1u << SUB_WARP_SIZE) - 1) << (sub_warp_id * SUB_WARP_SIZE);
    
    // Storage: LINEAR GAP only needs H (no E/F tracking)
    half2 row_H[TPT];
    
    // Boundaries: half2.x = seq0, half2.y = seq1 (LINEAR: only need H)
    half2 prev_boundary_H = __float2half2_rn(0.0f);
    half2 prev2_boundary_H = __float2half2_rn(0.0f);
    
    const half2 gap_penalty_h2 = __float2half2_rn((float)gap_penalty);
    const half2 ZERO_H2 = __float2half2_rn(0.0f);
    
    #pragma unroll
    for (int k = 0; k < TPT; ++k) row_H[k] = ZERO_H2;

    half2 max_score = ZERO_H2;  // Track max for both sequences
    
    // TRANSPOSED: Iterate over query positions instead of DB positions
    const int total_iters = query_seq_len + SUB_WARP_SIZE - 1;  // Use SUB_WARP_SIZE for wavefront

    // Load DB characters for each thread's TPT positions (TRANSPOSED!)
    // Precompute COMBINED index for 3D BLOSUM lookup from shared memory
    // Layout: [db0][db1][qchar], so we precompute db0*N*N + db1*N
    int db_combined_base[TPT];
    half2 mask[TPT];  // Validity mask: 1.0 if in-bounds, 0.0 if out-of-bounds
    
    #pragma unroll
    for (int k = 0; k < TPT; ++k) {
        int i_db = db_begin + k;

        int8_t c0, c1;
        // Use padding character (index 20) for out-of-bounds positions
        if (i_db < length0) {
            c0 = char_to_uint[ascii0[i_db]];
            if (c0 < 0) c0 = 20; // Clamp negative indices to padding 
        } else {
            c0 = 20;  // Padding character
        }
        if (i_db < length1) {
            c1 = char_to_uint[ascii1[i_db]];
            if (c1 < 0) c1 = 20; // Clamp negative indices to padding
        } else {
            c1 = 20;  // Padding character
        }
        
        // Precompute base index: db0 * N*N + db1 * N (then add qchar in the loop)
        db_combined_base[k] = c0 * NUM_AMINO_ACIDS_CUDA_2 + c1 * NUM_AMINO_ACIDS_CUDA;
    }

    // No preloading - will read directly from shared memory in the loop
    
    // Partial unroll of outer loop for better instruction scheduling
    #pragma unroll 16
    for (int iter = 0; iter < total_iters; ++iter) {
        int j_query = iter - sub_lane;  // Use sub_lane for sub-warp processing
        bool active_query = (j_query >= 0 && j_query < query_seq_len);
        
        // Check if we have any active database positions for this sub-warp
        int max_db_pos = db_begin + TPT - 1;
        bool active_db = (db_begin < length0 || db_begin < length1);  // At least one sequence has data in this range
        
        // Shuffle boundaries (half2 contains both sequences) - use sub-warp mask
        // LINEAR GAP: only need H boundaries
        half2 diag_H = __shfl_up_sync(sub_warp_mask, prev2_boundary_H, 1, SUB_WARP_SIZE);
        half2 left_boundary_H = __shfl_up_sync(sub_warp_mask, prev_boundary_H, 1, SUB_WARP_SIZE);
        
        // First thread in sub-warp gets zero
        if (sub_lane == 0) {
            diag_H = ZERO_H2;
            left_boundary_H = ZERO_H2;
        }

        // Skip iteration if query position invalid OR no active database positions
        if (!active_query || !active_db) continue;
        
        if (lane == 0) {
            diag_H = ZERO_H2;  // LINEAR: reset for lane 0
            left_boundary_H = ZERO_H2;
        }

        // All threads get the SAME query character for this iteration
        int qchar = query_seq_indices[j_query];

        // BLOSUM lookup from shared memory
        const half2* s_blosum_ptr = s_blosum_combined + qchar;
        
        // LINEAR GAP: Only track H values (no separate E/F tracking)
        half2 h_left_H = left_boundary_H;
        half2 prev_h_up = diag_H;  // Initialize to diag_H - eliminates conditional!

        // Process all TPT positions
        #pragma unroll
        for (int k = 0; k < TPT; ++k) {
            // BLOSUM lookup from shared memory
            half2 sub_packed = s_blosum_ptr[db_combined_base[k]];
            
            // Diagonal H for M computation (no conditional needed!)
            half2 h_diag_H = prev_h_up;
            
            // Load H from previous row (for up gap)
            half2 h_up_H = row_H[k];
            
            // LINEAR GAP: H = max(0, diag + sub, left_H - gap, up_H - gap)
            half2 M_term = __hadd2(h_diag_H, sub_packed);  // 1 add (match/mismatch)
            half2 left_gap = __hsub2(h_left_H, gap_penalty_h2);  // 1 sub (horizontal gap)
            half2 up_gap = __hsub2(h_up_H, gap_penalty_h2);  // 1 sub (vertical gap)
            
            half2 cur_H = __hmax2(__hmax2(__hmax2(
                M_term,           // Match/mismatch
                left_gap),        // Horizontal gap (from left)
                up_gap),          // Vertical gap (from up)
                ZERO_H2);         // Local alignment: max with 0
            
            // Track maximum score
            max_score = __hmax2(max_score, cur_H);
            
            // Save h_up for diagonal in next k
            prev_h_up = row_H[k];
            
            // Store for next iteration
            row_H[k] = cur_H;
        
            // Update left value for next k iteration (within same row)
            h_left_H = cur_H;
        }
                
        // Update boundaries for next row (LINEAR: only need H)
        prev2_boundary_H = prev_boundary_H;
        prev_boundary_H = row_H[TPT-1];
    }

    // Reduce within sub-warp for both sequences
    __shared__ half2 s_reduction_array[BLOCK_DIM_X];
    s_reduction_array[threadIdx.x] = max_score;
    __syncthreads();

    // Each sub-warp reduces independently within its SUB_WARP_SIZE threads
    for (int offset = SUB_WARP_SIZE / 2; offset > 0; offset >>= 1) {
        if (sub_lane < offset) {
            s_reduction_array[threadIdx.x] = __hmax2(s_reduction_array[threadIdx.x], s_reduction_array[threadIdx.x + offset]);
        }
        __syncthreads();
    }

    // Each sub-warp writes its own results (sub_lane 0 of each sub-warp)
    if (sub_lane == 0) {
        half2 final_max = s_reduction_array[threadIdx.x];
        if (seq_valid) {
            out_scores[seq0_idx] = __half2int_rn(__low2half(final_max));
        }
        if (seq1_valid && seq1_idx < num_db_seqs) {
            out_scores[seq1_idx] = __half2int_rn(__high2half(final_max));
        }
    }

#endif
}

