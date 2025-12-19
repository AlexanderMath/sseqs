/*
    Profile-based Smith-Waterman using SIMD half2.
    
    Uses a Position-Specific Scoring Matrix (PSSM) instead of BLOSUM62.
    PSSM shape: (query_len x 20) where each row is position-specific scores.
    
    Layout: pssm[query_position * 20 + amino_acid_index]
*/
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdint.h>
#include <cuda_fp16.h>

// Import COMPILE_PROFILE flag and char_to_uint from sw.cu
#ifndef COMPILE_PROFILE
#define COMPILE_PROFILE 1
#endif

// char_to_uint is defined in sw.cu - declare as extern
//extern __device__ __constant__ int8_t char_to_uint[256];

template<int BLOCK_DIM_X, int SUB_WARP_SIZE, int TPT>
__global__ void sw_kernel_profile(
    const int8_t* pssm,              // PSSM matrix: (query_len x 20), int8_t values
    int query_seq_len,
    int* good_idx,
    const uint8_t* ascii,
    int* starts,
    int num_db_seqs,
    int* out_scores,
    int gap_open,
    int gap_extend) {
    
#if COMPILE_PROFILE
    
    // DOUBLE-PACKING with SUB-WARPS: Each sub-warp processes 2 database sequences
    // half2.x = sequence 0, half2.y = sequence 1 (different sequences, same query position!)
    constexpr int WARP_SIZE = 32;
    constexpr int NUM_AA = 20;  // Standard amino acids
    
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
    
    // TRANSPOSED: Iterate over query positions instead of DB positions
    const int total_iters = query_seq_len + SUB_WARP_SIZE - 1;  // Use SUB_WARP_SIZE for wavefront

    // Load DB characters for each thread's TPT positions (TRANSPOSED!)
    // Store target amino acid indices for PSSM lookup
    int8_t db_char0[TPT];  // Target AA index for seq0
    int8_t db_char1[TPT];  // Target AA index for seq1
    
    #pragma unroll
    for (int k = 0; k < TPT; ++k) {
        int i_db = db_begin + k;

        // Use index 20 (maps to 0 score) for out-of-bounds positions
        if (i_db < length0) {
            db_char0[k] = char_to_uint[ascii0[i_db]];
            if (db_char0[k] < 0 || db_char0[k] >= NUM_AA) db_char0[k] = 0; // Clamp invalid to Ala
        } else {
            db_char0[k] = 0;  // Out of bounds -> score 0 anyway
        }
        if (i_db < length1) {
            db_char1[k] = char_to_uint[ascii1[i_db]];
            if (db_char1[k] < 0 || db_char1[k] >= NUM_AA) db_char1[k] = 0; // Clamp invalid to Ala
        } else {
            db_char1[k] = 0;  // Out of bounds -> score 0 anyway
        }
    }
    
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

        // All threads use the SAME query POSITION for this iteration
        // PSSM row pointer for this query position
        const int8_t* pssm_row = pssm + j_query * NUM_AA;
        
        // E = horizontal gap (rolling scalar), F = vertical gap (array)
        half2 h_left_H = left_boundary_H;
        half2 h_left_E = left_boundary_E;  // E is a rolling scalar
        //half2 prev_h_up = ZERO_H2;  // For tracking diagonal values across k iterations
        half2 prev_h_up = diag_H;  // For tracking diagonal values across k iterations

        // Check if all TPT positions are in bounds for both sequences
        #pragma unroll
        for (int k = 0; k < TPT; ++k) {
            // PSSM lookup: score = pssm[query_pos][target_aa]
            // Pack scores for both sequences into half2
            int8_t score0 = pssm_row[db_char0[k]];
            int8_t score1 = pssm_row[db_char1[k]];
            half2 sub_packed = __halves2half2(__int2half_rn(score0), __int2half_rn(score1));
            
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
            
            // NEW: Update E and F using CURRENT H
            // E = max(E - gap_extend, H - gap_open)
            half2 cur_E = __hmax2(__hsub2(h_left_E, gap_extend_h2), H_minus_open);  // 1 sub + 1 max
            
            // F = max(F - gap_extend, H - gap_open)
            half2 cur_F = __hmax2(__hsub2(h_up_F, gap_extend_h2), H_minus_open);  // 1 sub + 1 max
            
            // Track maximum only for valid positions (mask out-of-bounds)
            /*int i_db = db_begin + k;
            half masked_H_x = (i_db < length0) ? __low2half(cur_H) : __float2half(-10000.0f);
            half masked_H_y = (i_db < length1) ? __high2half(cur_H) : __float2half(-10000.0f);
            half2 masked_H = __halves2half2(masked_H_x, masked_H_y);
            max_score = __hmax2(max_score, masked_H);*/
            max_score = __hmax2(max_score, cur_H);
            
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