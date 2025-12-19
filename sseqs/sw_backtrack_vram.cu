
// Direction encoding for traceback (2 bits per sequence)
#define DIR_ZERO   0  // H came from zero (start of local alignment)
#define DIR_MATCH  1  // H came from diagonal (match/mismatch)  
#define DIR_INSERT 2  // H came from E (insertion/left)
#define DIR_DELETE 3  // H came from F (deletion/up)

template<int BLOCK_DIM_X, int SUB_WARP_SIZE, int TPT>
__global__ void sw_kernel_affine_backtrack(
    const unsigned char* query_seq_indices,
    int query_seq_len,
    int* good_idx,
    const uint8_t* ascii,
    int* starts,
    int num_db_seqs,
    int* out_scores,
    int* out_alignment_lens,  // Output: length of each alignment
    int* out_end_i,           // Output: ending i position (target)
    int* out_end_j,           // Output: ending j position (query)
    char* out_alignment_ops,  // Output: alignment operations (M/I/D), shape [num_db_seqs, max_len]
    int max_align_len,        // Maximum alignment length per sequence
    int gap_open,
    int gap_extend,
    uint8_t* global_directions, // Global memory: [num_seqs][m][n] direction pointers (4 bits per cell)
    int max_target_len,       // Maximum target sequence length (m)
    int max_query_len         // Maximum query sequence length (n)
) {
    
#if COMPILE_BACKTRACK
    

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
    
    // Load 3D combined BLOSUM into shared memory FIRST (before any early returns!)
    // Layout: [db0][db1][qchar] for better access pattern
    __shared__ half2 s_blosum_combined[NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA];
    
    // Cooperative load - ALL threads must participate before __syncthreads()
    for (int i = threadIdx.x; i < NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA; i += BLOCK_DIM_X) {
        s_blosum_combined[i] = blosum62_combined_cuda_global[i];
    }
    __syncthreads();
    
    // NOW we can do bounds checks (NO early returns - more __syncthreads() later!)
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

    // TODO REMOVE:THIS MAKES KERNEL 30-40% SLOWER. 
    #if PRINT_DEBUG
    if (length0 != length1) printf("length0 != length1: %d %d\n", length0, length1);
    #endif

    // Mark invalid if this thread's entire range is out of bounds (NO early return!)
    if (!seq_valid || (db_begin >= length0 && db_begin >= length1)) {
        length0 = 0;
        length1 = 0;
    }

    // Compute sub-warp mask for shuffle operations (only communicate within sub-warp)
    const unsigned int sub_warp_mask = ((1u << SUB_WARP_SIZE) - 1) << (sub_warp_id * SUB_WARP_SIZE);
    
    // Storage: H and F(vertical gap) as arrays; E(horizontal gap) is a rolling scalar
    half2 row_H[TPT], row_F[TPT];
    
    // Boundaries: half2.x = seq0, half2.y = seq1
    half2 prev_boundary_H = __float2half2_rn(0.0f);
    half2 prev_boundary_E = __float2half2_rn(0.0f);
    half2 prev_boundary_F = __float2half2_rn(0.0f);
    half2 prev2_boundary_H = __float2half2_rn(0.0f);
    half2 prev2_boundary_E = __float2half2_rn(0.0f);
    half2 prev2_boundary_F = __float2half2_rn(0.0f);

    half2 z = __float2half2_rn(10.0f); // todo: debug variable, remove in final version when all is done. 
    
    const half2 gap_open_h2 = __float2half2_rn((float)gap_open);
    const half2 gap_extend_h2 = __float2half2_rn((float)gap_extend);
    const half2 ZERO_H2 = __float2half2_rn(0.0f);
    
    //for (int k = 0; k < 0; ++k) { // 2.01 -> 2.15 TCUPs
    #pragma unroll
    for (int k = 0; k < TPT; ++k) row_H[k] = row_F[k] = ZERO_H2;

    half2 max_score = ZERO_H2;  // Track max for both sequences
    int max_i_0 = 0, max_j_0 = 0;  // Position of max score for seq0
    int max_i_1 = 0, max_j_1 = 0;  // Position of max score for seq1
    
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
        half2 diag_H = __shfl_up_sync(sub_warp_mask, prev2_boundary_H, 1, SUB_WARP_SIZE);
        half2 diag_E = __shfl_up_sync(sub_warp_mask, prev2_boundary_E, 1, SUB_WARP_SIZE);
        half2 diag_F = __shfl_up_sync(sub_warp_mask, prev2_boundary_F, 1, SUB_WARP_SIZE);
        half2 left_boundary_H =  __shfl_up_sync(sub_warp_mask, prev_boundary_H, 1, SUB_WARP_SIZE);
        half2 left_boundary_E = __shfl_up_sync(sub_warp_mask, prev_boundary_E, 1, SUB_WARP_SIZE);
        half2 left_boundary_F = __shfl_up_sync(sub_warp_mask, prev_boundary_F, 1, SUB_WARP_SIZE);
        
        // First thread in sub-warp gets zero
        if (sub_lane == 0) {
            diag_H = diag_E = diag_F = ZERO_H2;
            left_boundary_H = left_boundary_E = left_boundary_F = ZERO_H2;
        }

        // Skip iteration if query position invalid OR no active database positions
        if (!active_query || !active_db) continue;
        
        if (lane == 0) {
            diag_H = diag_E = diag_F = ZERO_H2;
            left_boundary_H = left_boundary_E = left_boundary_F = ZERO_H2;
        }

        // All threads get the SAME query character for this iteration
        int qchar = query_seq_indices[j_query];

        // BLOSUM lookup from shared memory - precompute base pointer to avoid repeated adds
        const half2* s_blosum_ptr = s_blosum_combined + qchar;
        
        // E = horizontal gap (rolling scalar), F = vertical gap (array)
        half2 h_left_H = left_boundary_H;
        half2 h_left_E = left_boundary_E;  // E is a rolling scalar
        //half2 prev_h_up = ZERO_H2;  // For tracking diagonal values across k iterations
        half2 prev_h_up = diag_H;  // For tracking diagonal values across k iterations

        // Check if all TPT positions are in bounds for both sequences
        #pragma unroll
        for (int k = 0; k < TPT; ++k) {
            // BLOSUM lookup from shared memory (using precomputed pointer)
            // Note: db_combined_base uses padding char (index 20) for out-of-bounds, which scores -1
            half2 sub_packed = s_blosum_ptr[db_combined_base[k]];
            
            // Load from previous row
            half2 h_up_F = row_F[k];  // F from previous row
            
            // Diagonal H for M computation
            //half2 h_diag_H = (k == 0) ? diag_H : prev_h_up; // <- gemmini thinks compiler is clever enough to remove with #pragma unroll
            half2 h_diag_H = prev_h_up; // <- gemmini thinks compiler is clever enough to remove with #pragma unroll
            
            // NEW: Compute H FIRST using previous E and F
            // H = max(0, diag + sub, E, F)
            half2 M_term = __hadd2(h_diag_H, sub_packed);  // 1 add
            half2 cur_H = __hmax2(__hmax2(__hmax2(
                M_term,           // Match/mismatch
                h_left_E),        // Horizontal gap (left)
                h_up_F),          // Vertical gap (up)
                ZERO_H2);         // 3 max
            
            // NEW: Precompute H - gap_open (KEY OPTIMIZATION!)
            half2 H_minus_open = __hsub2(cur_H, gap_open_h2);  // 1 sub
            
            // ============================================================================
            // STORE DIRECTION POINTERS TO GLOBAL MEMORY (for backtracking later)
            // ============================================================================
            // Storage: 1 byte per cell, containing 2x2-bit direction pointers
            // Lower 2 bits: direction for seq0 (.x), Upper 2 bits: direction for seq1 (.y)
            // IMPORTANT: Compare H with the OLD E/F values that were used to compute H!
            
            int i_db = db_begin + k;  // Target (database) position
            int seq_pair_idx = seq0_idx;  // Base index for this pair
            
            // Determine direction for seq0 (.x component)
            // Direction indicates how we ARRIVED at this cell
            half H_x = __low2half(cur_H);
            half E_x = __low2half(h_left_E);  // OLD E value (before update)
            half F_x = __low2half(h_up_F);    // OLD F value (before update)
            half M_x = __low2half(M_term);
            half zero = __float2half(0.0f);
            
            uint8_t dir0;
            // Check in order of preference: zero, match, insert, delete
            if (__hle(H_x, zero)) {  // H <= 0 means came from zero
                dir0 = DIR_ZERO;
            } else if (__heq(H_x, M_x)) {  // Check if came from match
                dir0 = DIR_MATCH;
            } else if (__heq(H_x, E_x)) {  // Check if came from insertion
                dir0 = DIR_INSERT;
            } else {  // Must have came from deletion
                dir0 = DIR_DELETE;
            }
            
            // Determine direction for seq1 (.y component)
            half H_y = __high2half(cur_H);
            half E_y = __high2half(h_left_E);  // OLD E value
            half F_y = __high2half(h_up_F);    // OLD F value
            half M_y = __high2half(M_term);
            
            uint8_t dir1;
            if (__hle(H_y, zero)) {
                dir1 = DIR_ZERO;
            } else if (__heq(H_y, M_y)) {
                dir1 = DIR_MATCH;
            } else if (__heq(H_y, E_y)) {
                dir1 = DIR_INSERT;
            } else {
                dir1 = DIR_DELETE;
            }
            
            // Pack both directions into 1 byte (4 bits each, only using lower 2 bits each)
            uint8_t packed_dir = (dir1 << 2) | dir0;
            
            // Write to global memory
            // Note: seq_pair_idx is the base seq index (0, 2, 4...), convert to pair index
            int pair_idx = seq_pair_idx / 2;
            size_t global_idx = (size_t)pair_idx * max_query_len * max_target_len + 
                               j_query * max_target_len + i_db;
            global_directions[global_idx] = packed_dir;
            
            // DEBUG: Print for first pair only, at max score position
            #if PRINT_DEBUG
            if (pair_idx == 0 && j_query < 3 && i_db < 3) {
                printf("FWD pair=%d seq0=%d seq1=%d i=%d j=%d | H0=%.1f H1=%.1f | dir0=%d dir1=%d | idx=%llu packed=0x%02x\n",
                       pair_idx, seq0_idx, seq1_idx, i_db, j_query,
                       __half2float(H_x), __half2float(H_y),
                       dir0, dir1, (unsigned long long)global_idx, packed_dir);
            }
            #endif
            // ============================================================================
            
            // NOW update E and F for the NEXT iteration
            // E = max(E - gap_extend, H - gap_open)
            half2 cur_E = __hmax2(__hsub2(h_left_E, gap_extend_h2), H_minus_open);  // 1 sub + 1 max
            
            // F = max(F - gap_extend, H - gap_open)
            half2 cur_F = __hmax2(__hsub2(h_up_F, gap_extend_h2), H_minus_open);  // 1 sub + 1 max
            
            // Track maximum and its position
            H_x = __low2half(cur_H);
            H_y = __high2half(cur_H);
            half max_x = __low2half(max_score);
            half max_y = __high2half(max_score);
            
            // Update max for seq0
            if (__hgt(H_x, max_x)) {
                max_x = H_x;
                max_i_0 = i_db;
                max_j_0 = j_query;
            }
            
            // Update max for seq1
            if (__hgt(H_y, max_y)) {
                max_y = H_y;
                max_i_1 = i_db;
                max_j_1 = j_query;
            }
            
            max_score = __halves2half2(max_x, max_y);
            
            // Save h_up for diagonal in next k
            prev_h_up = row_H[k];
            
            // Store for next iteration
            row_H[k] = cur_H; 
            row_F[k] = cur_F;
        
            // Update left values for next k iteration (within same row)
            h_left_H = cur_H;
            h_left_E = cur_E;  // E is a rolling scalar

        }
                
        // Update boundaries for next row
        prev2_boundary_H = prev_boundary_H;
        prev2_boundary_E = prev_boundary_E;
        prev2_boundary_F = prev_boundary_F;
        prev_boundary_H = row_H[TPT-1];
        prev_boundary_E = h_left_E;  // Use the final scalar E value
        prev_boundary_F = row_F[TPT-1];
    }

    // Reduce within sub-warp for both sequences
    __shared__ half2 s_reduction_array[BLOCK_DIM_X];
    __shared__ int s_max_i_0[BLOCK_DIM_X], s_max_j_0[BLOCK_DIM_X];
    __shared__ int s_max_i_1[BLOCK_DIM_X], s_max_j_1[BLOCK_DIM_X];
    
    s_reduction_array[threadIdx.x] = max_score;
    s_max_i_0[threadIdx.x] = max_i_0;
    s_max_j_0[threadIdx.x] = max_j_0;
    s_max_i_1[threadIdx.x] = max_i_1;
    s_max_j_1[threadIdx.x] = max_j_1;
    __syncthreads();

    // Each sub-warp reduces independently within its SUB_WARP_SIZE threads
    for (int offset = SUB_WARP_SIZE / 2; offset > 0; offset >>= 1) {
        if (sub_lane < offset) {
            half2 my_max = s_reduction_array[threadIdx.x];
            half2 other_max = s_reduction_array[threadIdx.x + offset];
            
            // Update seq0 max
            if (__hgt(__low2half(other_max), __low2half(my_max))) {
                my_max = __halves2half2(__low2half(other_max), __high2half(my_max));
                s_max_i_0[threadIdx.x] = s_max_i_0[threadIdx.x + offset];
                s_max_j_0[threadIdx.x] = s_max_j_0[threadIdx.x + offset];
            }
            
            // Update seq1 max
            if (__hgt(__high2half(other_max), __high2half(my_max))) {
                my_max = __halves2half2(__low2half(my_max), __high2half(other_max));
                s_max_i_1[threadIdx.x] = s_max_i_1[threadIdx.x + offset];
                s_max_j_1[threadIdx.x] = s_max_j_1[threadIdx.x + offset];
            }
            
            s_reduction_array[threadIdx.x] = my_max;
        }
        __syncthreads();
    }

    // Each sub-warp writes its own results (sub_lane 0 of each sub-warp)
    if (sub_lane == 0) {
        half2 final_max = s_reduction_array[threadIdx.x];
        if (seq_valid) {
            out_scores[seq0_idx] = __half2int_rn(__low2half(final_max));
            out_end_i[seq0_idx] = s_max_i_0[threadIdx.x];
            out_end_j[seq0_idx] = s_max_j_0[threadIdx.x];
            
            #if PRINT_DEBUG
            if (seq0_idx < 2) {
                printf("FWD_RESULT seq0=%d score=%d end_i=%d end_j=%d\n",
                       seq0_idx, __half2int_rn(__low2half(final_max)),
                       s_max_i_0[threadIdx.x], s_max_j_0[threadIdx.x]);
            }
            #endif
        }
        if (seq1_valid && seq1_idx < num_db_seqs) {
            out_scores[seq1_idx] = __half2int_rn(__high2half(final_max));
            out_end_i[seq1_idx] = s_max_i_1[threadIdx.x];
            out_end_j[seq1_idx] = s_max_j_1[threadIdx.x];
            
            #if PRINT_DEBUG
            if (seq1_idx < 2) {
                printf("FWD_RESULT seq1=%d score=%d end_i=%d end_j=%d\n",
                       seq1_idx, __half2int_rn(__high2half(final_max)),
                       s_max_i_1[threadIdx.x], s_max_j_1[threadIdx.x]);
            }
            #endif
        }
    }

    // Set alignment lengths to 0 (will be filled by traceback kernel)
    if (sub_lane == 0) {
        if (seq_valid) {
            out_alignment_lens[seq0_idx] = 0;
        }
        if (seq1_valid) {
            out_alignment_lens[seq1_idx] = 0;
        }
    }


#endif
}

// ============================================================================
// TRACEBACK KERNEL: Follow direction pointers to reconstruct alignment
// ============================================================================
// Simple kernel: 1 thread per sequence, serial traceback
__global__ void sw_traceback_kernel(
    const uint8_t* global_directions,  // Direction pointers from forward pass
    const int* scores,                  // Scores to find max position
    int* out_alignment_lens,            // Output: length of each alignment
    int* out_end_i,                     // Output: ending i position (target)
    int* out_end_j,                     // Output: ending j position (query)
    char* out_alignment_ops,            // Output: alignment operations (M/I/D)
    int max_align_len,                  // Maximum alignment length
    int num_seqs,
    int max_target_len,
    int max_query_len
) {
    int seq_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (seq_idx >= num_seqs) return;
    
    // Determine if this is seq0 or seq1 in a pair
    // Forward pass stores with seq_pair_idx = seq0_idx (0, 2, 4, 6...)
    // So seq_idx=0,1 are in pair 0; seq_idx=2,3 in pair 2, etc.
    int seq_pair_idx = (seq_idx / 2) * 2;  // Round down to even number
    bool is_seq1 = (seq_idx % 2) == 1;     // Odd indices are seq1
    
    // Get max score position from forward pass outputs
    int i = out_end_i[seq_idx];
    int j = out_end_j[seq_idx];
    
    #if PRINT_DEBUG
    if (seq_idx < 2) {
        int pair_idx = seq_pair_idx / 2;
        printf("TRACE_START seq=%d pair_idx=%d seq_pair_idx=%d is_seq1=%d start_i=%d start_j=%d max_target=%d max_query=%d\n",
               seq_idx, pair_idx, seq_pair_idx, is_seq1, i, j, max_target_len, max_query_len);
    }
    #endif
    
    // Traceback: follow direction pointers
    int align_len = 0;
    char* alignment = out_alignment_ops + seq_idx * max_align_len;
    
    while (i > 0 && j > 0 && align_len < max_align_len) {
        // Read direction pointer using pair_idx
        // Note: seq_pair_idx is the base seq index (0, 2, 4...), convert to pair index
        int pair_idx = seq_pair_idx / 2;
        size_t global_idx = (size_t)pair_idx * max_query_len * max_target_len + 
                           j * max_target_len + i;
        uint8_t packed_dir = global_directions[global_idx];
        
        // Extract direction: lower 2 bits for seq0, upper 2 bits for seq1
        uint8_t dir = is_seq1 ? ((packed_dir >> 2) & 0x3) : (packed_dir & 0x3);
        
        #if PRINT_DEBUG
        if (seq_idx < 2 && align_len < 10) {
            printf("TRACE seq=%d step=%d i=%d j=%d | idx=%llu packed=0x%02x dir=%d | ",
                   seq_idx, align_len, i, j, (unsigned long long)global_idx, packed_dir, dir);
        }
        #endif
        
        if (dir == DIR_ZERO) {
            #if PRINT_DEBUG
            if (seq_idx < 2 && align_len < 10) printf("ZERO (stop)\n");
            #endif
            break;  // Reached start of local alignment
        } else if (dir == DIR_MATCH) {
            alignment[align_len++] = 'M';
            #if PRINT_DEBUG
            if (seq_idx < 2 && align_len <= 10) printf("MATCH\n");
            #endif
            i--;
            j--;
        } else if (dir == DIR_INSERT) {
            alignment[align_len++] = 'I';
            #if PRINT_DEBUG
            if (seq_idx < 2 && align_len <= 10) printf("INSERT\n");
            #endif
            j--;  // Move left (insertion in query)
        } else if (dir == DIR_DELETE) {
            alignment[align_len++] = 'D';
            #if PRINT_DEBUG
            if (seq_idx < 2 && align_len <= 10) printf("DELETE\n");
            #endif
            i--;  // Move up (deletion in query)
        }
    }
    
    // Store results
    // IMPORTANT: Keep the ending positions (i_end, j_end) unchanged!
    // They were set by forward pass and Python needs them for reconstruction
    out_alignment_lens[seq_idx] = align_len;
    // Do NOT overwrite out_end_i and out_end_j - Python needs the original ending positions!
    // out_end_i[seq_idx] = i;  // This would overwrite ending with starting position
    // out_end_j[seq_idx] = j;  // This would overwrite ending with starting position
    
    #if PRINT_DEBUG
    if (seq_idx < 2) {
        printf("TRACE_END seq=%d align_len=%d start_i=%d start_j=%d (end_i=%d end_j=%d unchanged)\n",
               seq_idx, align_len, i, j, out_end_i[seq_idx], out_end_j[seq_idx]);
    }
    #endif
    
    // Reverse alignment (we traced backwards)
    for (int k = 0; k < align_len / 2; k++) {
        char tmp = alignment[k];
        alignment[k] = alignment[align_len - 1 - k];
        alignment[align_len - 1 - k] = tmp;
    }
}