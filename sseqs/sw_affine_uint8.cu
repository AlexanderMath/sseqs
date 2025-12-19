template<int BLOCK_DIM_X, int SUB_WARP_SIZE, int TPT>
__global__ void sw_kernel_affine_int8( const unsigned char* query_seq_indices, int query_seq_len, int* good_idx, const uint8_t* ascii,
    int* starts, // call length earlier -- refactor to call starts, it contains the starting position of all sequences!
    int num_db_seqs, int* out_scores, int gap_open, int gap_extend) {
    
#if COMPILE_AFFINE_INT8
    
    // QUAD-PACKING with SUB-WARPS: Each sub-warp processes 4 database sequences using __vadd4
    // unsigned int packs 4× uint8 values = 4 sequences per thread
    constexpr int WARP_SIZE = 32;
    
    const int warp_id = threadIdx.x / WARP_SIZE;
    const int lane = threadIdx.x % WARP_SIZE;
    const int sub_warp_id = lane / SUB_WARP_SIZE;  // Which sub-warp within the warp
    const int sub_lane = lane % SUB_WARP_SIZE;      // Position within sub-warp
    
    constexpr int WARPS_PER_BLOCK = BLOCK_DIM_X / WARP_SIZE;
    constexpr int SUB_WARPS_PER_WARP = WARP_SIZE / SUB_WARP_SIZE;
    constexpr int SEQS_PER_SUB_WARP = 4;  // Quad-packing: 4 sequences per sub-warp
    constexpr int SEQS_PER_WARP = SUB_WARPS_PER_WARP * SEQS_PER_SUB_WARP;
    constexpr int SEQS_PER_BLOCK = WARPS_PER_BLOCK * SEQS_PER_WARP;
    
    // Each sub-warp processes 4 sequences
    const int seq0_idx = blockIdx.x * SEQS_PER_BLOCK + warp_id * SEQS_PER_WARP + sub_warp_id * 4;
    const int seq1_idx = seq0_idx + 1;
    const int seq2_idx = seq0_idx + 2;
    const int seq3_idx = seq0_idx + 3;
    
    
    // Get sequence info for all 4 sequences
    int start0 = starts[good_idx[seq0_idx]-1] + 1;
    int stop0 = starts[good_idx[seq0_idx]] + 1;
    int length0 = stop0 - start0 - 1;
    const uint8_t* ascii0 = ascii + start0;
    
    int start1 = starts[good_idx[seq1_idx]-1] + 1;
    int stop1 = starts[good_idx[seq1_idx]] + 1;
    int length1 = stop1 - start1 - 1;
    const uint8_t* ascii1 = ascii + start1;
    
    int start2 = starts[good_idx[seq2_idx]-1] + 1;
    int stop2 = starts[good_idx[seq2_idx]] + 1;
    int length2 = stop2 - start2 - 1;
    const uint8_t* ascii2 = ascii + start2;
    
    int start3 = starts[good_idx[seq3_idx]-1] + 1;
    int stop3 = starts[good_idx[seq3_idx]] + 1;
    int length3 = stop3 - start3 - 1;
    const uint8_t* ascii3 = ascii + start3;
    
    // Load 3D int8 BLOSUM into shared memory
    __shared__ uint16_t s_blosum_int8[NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA];
    
    // Cooperative load
    for (int i = threadIdx.x; i < NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA; i += BLOCK_DIM_X) {
        s_blosum_int8[i] = blosum62_int8_combined_global[i];
    }
    __syncthreads();
    
    // TRANSPOSED: Each thread handles TPT consecutive DB positions
    const int db_begin = sub_lane * TPT;
    
    // Compute sub-warp mask for shuffle operations
    const unsigned int sub_warp_mask = ((1u << SUB_WARP_SIZE) - 1) << (sub_warp_id * SUB_WARP_SIZE);
    
    // Storage: H and F(vertical gap) as arrays; E(horizontal gap) is a rolling scalar
    unsigned int row_H[TPT], row_F[TPT];
    
    // Boundaries: int packs 4 signed int8 sequences
    const unsigned int ZERO_I4 = 0x00000000u;  // Zero in each byte
    unsigned int prev_boundary_H = ZERO_I4;
    unsigned int prev_boundary_E = ZERO_I4;
    unsigned int prev_boundary_F = ZERO_I4;
    unsigned int prev2_boundary_H = ZERO_I4;
    unsigned int prev2_boundary_E = ZERO_I4;
    unsigned int prev2_boundary_F = ZERO_I4;
    
    // Gap penalties as packed signed int8
    unsigned int gap_open_i4 = ((gap_open & 0xFF) << 24) | ((gap_open & 0xFF) << 16) | ((gap_open & 0xFF) << 8) | (gap_open & 0xFF);
    unsigned int gap_extend_i4 = ((gap_extend & 0xFF) << 24) | ((gap_extend & 0xFF) << 16) | ((gap_extend & 0xFF) << 8) | (gap_extend & 0xFF);
        
    #pragma unroll
    for (int k = 0; k < TPT; ++k) row_H[k] = row_F[k] = ZERO_I4;

    unsigned int max_score = ZERO_I4;  // Track max for all 4 sequences
    unsigned int overflow_flag = 0u;    // Track if any sequence overflowed (bit per sequence)
    const unsigned int OVERFLOW_THRESHOLD = 0x78787878u;  // 120 in each byte - conservative threshold
    
    // TRANSPOSED: Iterate over query positions instead of DB positions
    const int total_iters = query_seq_len + SUB_WARP_SIZE - 1;

    // Load DB characters for each thread's TPT positions (TRANSPOSED!)
    // We need two base indices - one for seq0+seq1, one for seq2+seq3
    int db_combined_base01[TPT];  // For sequences 0 and 1
    int db_combined_base23[TPT];  // For sequences 2 and 3
    
    #pragma unroll
    for (int k = 0; k < TPT; ++k) {
        int i_db = db_begin + k;
        int8_t c0 = char_to_uint[ascii0[i_db]];
        int8_t c1 = char_to_uint[ascii1[i_db]];
        int8_t c2 = char_to_uint[ascii2[i_db]];
        int8_t c3 = char_to_uint[ascii3[i_db]];
        
        // Clamp negative indices
        if (c0 < 0) c0 = 0;
        if (c1 < 0) c1 = 0;
        if (c2 < 0) c2 = 0;
        if (c3 < 0) c3 = 0;
        
        // Precompute indices for 3D BLOSUM lookup (same layout as half2 version)
        db_combined_base01[k] = c0 * NUM_AMINO_ACIDS_CUDA_2 + c1 * NUM_AMINO_ACIDS_CUDA;
        db_combined_base23[k] = c2 * NUM_AMINO_ACIDS_CUDA_2 + c3 * NUM_AMINO_ACIDS_CUDA;
    }

    // Partial unroll of outer loop for better instruction scheduling
    #pragma unroll 16
    for (int iter = 0; iter < total_iters; ++iter) {
        int j_query = iter - sub_lane;
        bool active_query = (j_query >= 0 && j_query < query_seq_len);
        
        // Shuffle boundaries (unsigned int contains 4 sequences)
        unsigned int diag_H = __shfl_up_sync(sub_warp_mask, prev2_boundary_H, 1, SUB_WARP_SIZE);
        unsigned int diag_E = __shfl_up_sync(sub_warp_mask, prev2_boundary_E, 1, SUB_WARP_SIZE);
        unsigned int diag_F = __shfl_up_sync(sub_warp_mask, prev2_boundary_F, 1, SUB_WARP_SIZE);
        unsigned int left_boundary_H =  __shfl_up_sync(sub_warp_mask, prev_boundary_H, 1, SUB_WARP_SIZE);
        unsigned int left_boundary_E = __shfl_up_sync(sub_warp_mask, prev_boundary_E, 1, SUB_WARP_SIZE);
        unsigned int left_boundary_F = __shfl_up_sync(sub_warp_mask, prev_boundary_F, 1, SUB_WARP_SIZE);
        
        // First thread in sub-warp gets zero
        if (sub_lane == 0) {
            diag_H = diag_E = diag_F = ZERO_I4;
            left_boundary_H = left_boundary_E = left_boundary_F = ZERO_I4;
        }

        // Skip iteration if query position invalid
        if (!active_query) continue;
        
        if (lane == 0) {
            diag_H = diag_E = diag_F = ZERO_I4;
            left_boundary_H = left_boundary_E = left_boundary_F = ZERO_I4;
        }

        // All threads get the SAME query character for this iteration
        int qchar = query_seq_indices[j_query];

        // BLOSUM lookup: load 2× uint16 (2x int8 each) and pack into uint (4x int8)
        // Layout is [db0][db1][qchar], so we need to compute offset differently
        
        // E = horizontal gap (rolling scalar), F = vertical gap (array)
        unsigned int h_left_H = left_boundary_H;
        unsigned int h_left_E = left_boundary_E;
        unsigned int prev_h_up = ZERO_I4;
        
        #pragma unroll
        for (int k = 0; k < TPT; ++k) {
            // Compute BLOSUM indices: layout is [db0][db1][qchar]
            uint16_t sub_packed01 = s_blosum_int8[db_combined_base01[k] + qchar];
            uint16_t sub_packed23 = s_blosum_int8[db_combined_base23[k] + qchar];
            
            // Extract int8 scores (already in correct format, just repack into uint)
            unsigned int sub_packed = ((sub_packed23 & 0xFFFF) << 16) | (sub_packed01 & 0xFFFF);
            
            // Load from previous row
            unsigned int h_up_F = row_F[k];
            
            // Diagonal H for M computation
            unsigned int h_diag_H = (k == 0) ? diag_H : prev_h_up;
            
            // Compute H using signed int8 SIMD operations
            // H = max(0, diag + sub, E, F)
            unsigned int M_term = __vadd4(h_diag_H, sub_packed);  // 4× int8 add
            //unsigned int M_term =  sub_packed;  // 4× int8 add
            unsigned int cur_H = __vmaxs4(__vmaxs4(__vmaxs4(
                M_term,           // Match/mismatch
                h_left_E),        // Horizontal gap (left)
                h_up_F),          // Vertical gap (up)
                ZERO_I4);         // 4× int8 max (clamps to 0)
            
            // Precompute H - gap_open
            unsigned int H_minus_open = __vsubss4(cur_H, gap_open_i4);  // 4× int8 sub (saturating)
            //unsigned int H_minus_open = cur_H;  // 4× int8 sub (saturating)
            
            // Update E and F using CURRENT H
            // E = max(E - gap_extend, H - gap_open)
            unsigned int cur_E = __vmaxs4(__vsubss4(h_left_E, gap_extend_i4), H_minus_open);
            //unsigned int cur_E = __vmaxs4(h_left_E, H_minus_open);
            //unsigned int cur_E = h_left_E;
            
            // F = max(F - gap_extend, H - gap_open)
            unsigned int cur_F = __vmaxs4(__vsubss4(h_up_F, gap_extend_i4), H_minus_open);
            //unsigned int cur_F = __vmaxs4(h_up_F, H_minus_open);
            //unsigned int cur_F = h_up_F;
            
            // Track maximum
            max_score = __vmaxs4(max_score, cur_H);
            //max_score =  cur_H;
            
            // Save h_up for diagonal in next k
            prev_h_up = row_H[k];
            
            // Store for next iteration
            row_H[k] = cur_H;
            row_F[k] = cur_F;
            
            // Update left values for next k iteration
            h_left_H = cur_H;
            h_left_E = cur_E;
        }
        
        // Update boundaries for next row
        prev2_boundary_H = prev_boundary_H;
        prev2_boundary_E = prev_boundary_E;
        prev2_boundary_F = prev_boundary_F;
        prev_boundary_H = row_H[TPT-1];
        prev_boundary_E = h_left_E;
        prev_boundary_F = row_F[TPT-1];
    }

    // Reduce within sub-warp for all 4 sequences
    __shared__ unsigned int s_reduction_array[BLOCK_DIM_X];
    __shared__ unsigned int s_overflow_flags[BLOCK_DIM_X];
    s_reduction_array[threadIdx.x] = max_score;
    s_overflow_flags[threadIdx.x] = overflow_flag;
    __syncthreads();

    // Each sub-warp reduces independently
    for (int offset = SUB_WARP_SIZE / 2; offset > 0; offset >>= 1) {
        if (sub_lane < offset) {
            s_reduction_array[threadIdx.x] = __vmaxs4(s_reduction_array[threadIdx.x], s_reduction_array[threadIdx.x + offset]);
            s_overflow_flags[threadIdx.x] |= s_overflow_flags[threadIdx.x + offset];
        }
        __syncthreads();
    }

    // Extract 4 scores from packed int (signed int8, no offset)
    // Mark overflows as -1 if overflow flag is set for that byte
    if (sub_lane == 0) {
        unsigned int final_max = s_reduction_array[threadIdx.x];
        unsigned int final_overflow = s_overflow_flags[threadIdx.x];
        
        int8_t score0 = (int8_t)((final_max      ) & 0xFF);
        int8_t score1 = (int8_t)((final_max >> 8 ) & 0xFF);
        int8_t score2 = (int8_t)((final_max >> 16) & 0xFF);
        int8_t score3 = (int8_t)((final_max >> 24) & 0xFF);
        
        // Check overflow flag for each byte
        bool overflow0 = (final_overflow & 0xFF) != 0;
        bool overflow1 = ((final_overflow >> 8) & 0xFF) != 0;
        bool overflow2 = ((final_overflow >> 16) & 0xFF) != 0;
        bool overflow3 = ((final_overflow >> 24) & 0xFF) != 0;
        
        out_scores[seq0_idx] = overflow0 ? -1 : score0;
        out_scores[seq1_idx] = overflow1 ? -1 : score1;
        out_scores[seq2_idx] = overflow2 ? -1 : score2;
        out_scores[seq3_idx] = overflow3 ? -1 : score3;
    }
#endif     
}