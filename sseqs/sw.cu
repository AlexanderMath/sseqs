#include <cuda.h>
#include <cuda_runtime.h>
#include <stdint.h> // For uint8_t
#include <cstdio>        

#define NUM_AMINO_ACIDS_CUDA 20
__device__ int8_t blosum62_matrix_cuda_global[NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA] = {
//   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
     4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0, // A
    -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3, // R
    -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, // N
    -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, // D
     0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1, // C
    -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, // Q
    -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, // E
     0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3, // G
    -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, // H
    -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3, // I
    -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1, // L
    -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, // K
    -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1, // M
    -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1, // F
    -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2, // P
     1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, // S
     0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0, // T
    -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3, // W
    -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1, // Y
     0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4  // V
};

__device__ __constant__ int8_t char_to_uint[256] = {
    [0 ... 255] = -1,  
    ['A'] =  0,
    ['R'] =  1,
    ['N'] =  2,
    ['D'] =  3,
    ['C'] =  4,
    ['Q'] =  5,
    ['E'] =  6,
    ['G'] =  7,
    ['H'] =  8,
    ['I'] =  9,
    ['L'] = 10,
    ['K'] = 11,
    ['M'] = 12,
    ['F'] = 13,
    ['P'] = 14,
    ['S'] = 15,
    ['T'] = 16,
    ['W'] = 17,
    ['Y'] = 18,
    ['V'] = 19
};

template<int BLOCK_DIM_X, int MAX_SHARED_QUERY_LEN_SW>
__global__ void sw_kernel(
    const unsigned char* query_seq_indices,
    int query_seq_len,
    int* good_idx,
    const uint8_t* ascii,
    int* starts, // array contains the starting position of sequences
    int num_db_seqs,
    int* out_scores,
    int gap_penalty
) {
    const int seq_idx_global = blockIdx.x;
    if (seq_idx_global >= num_db_seqs) return; 
    
    int start = starts[good_idx[seq_idx_global]-1] + 1;
    int stop = starts[good_idx[seq_idx_global]] + 1;
    int length = stop - start - 1;

    if (query_seq_len == 0 || length == 0 || query_seq_len > MAX_SHARED_QUERY_LEN_SW || stop < start) {
        //printf("sw_kernel skipping %d\n", seq_idx_global);
        out_scores[seq_idx_global] = -1;
        return; 
    }

    constexpr int TPT = MAX_SHARED_QUERY_LEN_SW / 32;
    ascii = ascii + start;  

    // Move BLOSUM62 matrix to shared memory 
    __shared__ int8_t s_blosum62_matrix[NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA]; // 400=0.4k
    for (int i = threadIdx.x; i < (NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA); i += BLOCK_DIM_X) {
        s_blosum62_matrix[i] = blosum62_matrix_cuda_global[i];
    }

    // Move query sequence to local variable. 
    __shared__ unsigned char s_query_seq_sdata[MAX_SHARED_QUERY_LEN_SW]; 
    for (int i = threadIdx.x; i < MAX_SHARED_QUERY_LEN_SW; i += BLOCK_DIM_X) {
        if (i < query_seq_len) s_query_seq_sdata[i] = query_seq_indices[i];
        else s_query_seq_sdata[i] = 0; 
    }

    const int lane   = threadIdx.x;      // 0-31 inside the war
    //const int warpId = threadIdx.x >> 5; // divide by 32 = 2**5 = num_threads. 
    int p_thread_max_score = 0; 

    /* ─────────────────────────  STRIPE  VERSION  ───────────────────────── */
    const int j_begin = lane * TPT;
    const int j_end   = min(j_begin + TPT, query_seq_len);   // open interval
    //uint16_t diag_val = 0; 
    uint16_t prev[TPT];// = {0};
    uint16_t curr[TPT];// = {0};
    #pragma unroll
    for (int k = 0; k < TPT; ++k) { prev[k] = 0; curr[k] = 0; }

    for (int i_db = 0; i_db < length; ++i_db)   {
        /* upper-left neighbour of the first column in my stripe */
        //int active_cols  = max(0,min(TPT, query_seq_len - j_begin));   // columns this lane owns
        uint16_t last_up = prev[TPT - 1];               // H(i-1, j_end-1)
        uint16_t diag_val = __shfl_up_sync(0xffffffff, last_up, 1);
        if (lane == 0) diag_val = 0;

        /* constant for the whole row */
        int8_t indx = ascii[i_db];
        int8_t db_idx = char_to_uint[indx]; // -1 when bad. 
        uint16_t h_left = 0;                    // H(i, j-1) inside my stripe

        /* ── column loop over my stripe ────────────────────────────────── */
        #pragma unroll
        for (int k = 0; k < TPT; ++k) {
            int j_query = j_begin + k;
            if (j_query >= query_seq_len) break;

            unsigned char query_char_val = s_query_seq_sdata[j_query]; // loading constant doesn't change speed here. 

            int8_t sub_score = (db_idx < 0 || db_idx >= NUM_AMINO_ACIDS_CUDA)
                                ? -4
                                :s_blosum62_matrix[db_idx * NUM_AMINO_ACIDS_CUDA + query_char_val]; // removing this makes 5.2 -> 4.69ms

            uint16_t h_up_pred   = prev[k];      // H(i-1, j)
            uint16_t h_diag_pred = diag_val;     // H(i-1, j-1)
            uint16_t h_left_pred = h_left;       // H(i,   j-1)

            // removing max to all adds makes it 5.2ms -> 3.96ms
            // rtx5090 has single operation for this stuff. 
            uint16_t cur = max(0,
                            max(h_diag_pred + sub_score,
                                        max(h_up_pred + gap_penalty,
                                                        h_left_pred + gap_penalty)
                                                        )
                                                        );
            curr[k] = cur;
            p_thread_max_score = max(p_thread_max_score, cur);

            diag_val = h_up_pred;   // becomes upper-left for next column
            h_left   = cur;         // carry within stripe
        }

        #pragma unroll
        for (int k = 0; k < TPT; ++k)  prev[k] = curr[k];   // next row's "up"
        diag_val = __shfl_up_sync(0xffffffff, prev[TPT-1], 1);
        if (lane == 0) diag_val = 0;
    }


    // Reduction to find max score in the block
    __shared__ int s_reduction_array[BLOCK_DIM_X]; 
    s_reduction_array[threadIdx.x] = p_thread_max_score;
    __syncwarp();

    // Parallel reduction in shared memory
    // Each iteration halves the number of active threads and values to compare
    for (int offset = BLOCK_DIM_X / 2; offset > 0; offset >>= 1) {
        // offset will be BLOCK_DIM_X/2, BLOCK_DIM_X/4, ..., 1
        if (threadIdx.x < offset) {
            s_reduction_array[threadIdx.x] = max(s_reduction_array[threadIdx.x], s_reduction_array[threadIdx.x + offset]);
        }
        __syncwarp(); 
    }

    // The final maximum score for the block is now in s_reduction_array[0]
    if (threadIdx.x == 0) {
        out_scores[seq_idx_global] = s_reduction_array[0];
    }
    __syncthreads();
    // 220ms
}

// penalizing opening. 
template<int BLOCK_DIM_X, int MAX_SHARED_QUERY_LEN_SW>
__global__ void sw_kernel_affine(
    const unsigned char* query_seq_indices,
    int query_seq_len,
    int* good_idx,
    const uint8_t* ascii,
    int* starts, // call length earlier -- refactor to call starts, it contains the starting position of all sequences!
    int num_db_seqs,
    int* out_scores,
    int gap_open,
    int gap_extend
) {
    const int seq_idx_global = blockIdx.x;
    if (seq_idx_global >= num_db_seqs) return;
    
    int start = starts[good_idx[seq_idx_global]-1] + 1;
    int stop = starts[good_idx[seq_idx_global]] + 1;
    int length = stop - start - 1;

    if (query_seq_len == 0 || length == 0 || query_seq_len > MAX_SHARED_QUERY_LEN_SW || stop < start) {
        out_scores[seq_idx_global] = -1;
        return;
    }

    constexpr int TPT = MAX_SHARED_QUERY_LEN_SW / 32;
    ascii = ascii + start;

    __shared__ int8_t s_blosum62_matrix[NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA];
    for (int i = threadIdx.x; i < (NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA); i += BLOCK_DIM_X) {
        s_blosum62_matrix[i] = blosum62_matrix_cuda_global[i];
    }

    __shared__ unsigned char s_query_seq_sdata[MAX_SHARED_QUERY_LEN_SW];
    for (int i = threadIdx.x; i < MAX_SHARED_QUERY_LEN_SW; i += BLOCK_DIM_X) {
        if (i < query_seq_len) s_query_seq_sdata[i] = query_seq_indices[i];
        else s_query_seq_sdata[i] = 0;
    }

    const int lane = threadIdx.x;
    //const int warpId = threadIdx.x >> 5;
    int p_thread_max_score = 0;

    const int j_begin = lane * TPT;
    const int j_end = min(j_begin + TPT, query_seq_len);
    // Change score type to int16_t for signed arithmetic
    const int16_t ZERO16 = 0;
    int16_t diag_val_M = 0, diag_val_Ix = 0, diag_val_Iy = 0;
    int16_t prev_M[TPT], prev_Ix[TPT], prev_Iy[TPT];
    int16_t curr_M[TPT], curr_Ix[TPT], curr_Iy[TPT];
    #pragma unroll
    for (int k = 0; k < TPT; ++k) {
        prev_M[k] = prev_Ix[k] = prev_Iy[k] = 0;
        curr_M[k] = curr_Ix[k] = curr_Iy[k] = 0;
    }

    for (int i_db = 0; i_db < length; ++i_db) {
        uint16_t last_up_M = prev_M[TPT - 1];
        uint16_t last_up_Ix = prev_Ix[TPT - 1];
        uint16_t last_up_Iy = prev_Iy[TPT - 1];
        diag_val_M = __shfl_up_sync(0xffffffff, last_up_M, 1);
        diag_val_Ix = __shfl_up_sync(0xffffffff, last_up_Ix, 1);
        diag_val_Iy = __shfl_up_sync(0xffffffff, last_up_Iy, 1);
        if (lane == 0) {
            diag_val_M = diag_val_Ix = diag_val_Iy = 0;
        }

        int8_t indx = ascii[i_db];
        int8_t db_idx = char_to_uint[indx];
        uint16_t h_left_M = 0, h_left_Ix = 0;
        //uint16_t h_left_Iy;// = 0;

        #pragma unroll
        for (int k = 0; k < TPT; ++k) {
            int j_query = j_begin + k;
            if (j_query >= query_seq_len) break;

            unsigned char query_char_val = s_query_seq_sdata[j_query];
            int8_t sub_score = (db_idx < 0 || db_idx >= NUM_AMINO_ACIDS_CUDA)
                                ? -4
                                : s_blosum62_matrix[db_idx * NUM_AMINO_ACIDS_CUDA + query_char_val];

            int16_t h_up_M  = prev_M[k];
            int16_t h_up_Ix = prev_Ix[k];
            int16_t h_up_Iy = prev_Iy[k];
            int16_t h_diag_M = diag_val_M;
            int16_t h_diag_Ix = diag_val_Ix;
            int16_t h_diag_Iy = diag_val_Iy;
            int16_t h_left_M_pred = h_left_M;
            int16_t h_left_Ix_pred = h_left_Ix;
            //int16_t h_left_Iy_pred = h_left_Iy;

            // M(i,j)
            int16_t cur_M = max(ZERO16,
                               max(h_diag_M + sub_score,
                                            max(h_diag_Ix + sub_score,
                                                         h_diag_Iy + sub_score)));
            // Ix(i,j) : gap in query (horizontal): coming from left
            int16_t cur_Ix = max(ZERO16,
                                max(h_left_M_pred - gap_open,
                                             h_left_Ix_pred - gap_extend));
            // Iy(i,j) : gap in target (vertical): coming from up
            int16_t cur_Iy = max(ZERO16,
                                max(h_up_M - gap_open,
                                             h_up_Iy - gap_extend));

            curr_M[k] = cur_M;
            curr_Ix[k] = cur_Ix;
            curr_Iy[k] = cur_Iy;

            p_thread_max_score = max(p_thread_max_score,
                                        max(cur_M, max(cur_Ix, cur_Iy)));

            diag_val_M = h_up_M;
            diag_val_Ix = h_up_Ix;
            diag_val_Iy = h_up_Iy;
            h_left_M = cur_M;
            h_left_Ix = cur_Ix;
            //h_left_Iy = cur_Iy;
        }

        #pragma unroll
        for (int k = 0; k < TPT; ++k) {
            prev_M[k] = curr_M[k];
            prev_Ix[k] = curr_Ix[k];
            prev_Iy[k] = curr_Iy[k];
        }
    }

    __shared__ int s_reduction_array[BLOCK_DIM_X];
    s_reduction_array[threadIdx.x] = p_thread_max_score;
    __syncwarp();

    for (int offset = BLOCK_DIM_X / 2; offset > 0; offset >>= 1) {
        if (threadIdx.x < offset) {
            s_reduction_array[threadIdx.x] = max(s_reduction_array[threadIdx.x], s_reduction_array[threadIdx.x + offset]);
        }
        __syncwarp();
    }

    if (threadIdx.x == 0) {
        out_scores[seq_idx_global] = s_reduction_array[0];
    }
    __syncthreads();
}



#define LAUNCH_CASE(LEN)                                                         \
    case LEN:                                                                    \
        sw_kernel<32, LEN><<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_val);                                                        \
        break;



extern "C"
void launch_sw_cuda_blosum62(
        const uint8_t* query_seq_indices_ptr,
        int            query_length,
        int*           good_idx,
        const uint8_t* ascii,
        int*           lengths,
        int            num_sequences_in_db, // Actual number of DB sequences
        int*           output_scores_ptr,
        int            gap_val,
        cudaStream_t   stream)
{
    dim3 num_blocks(num_sequences_in_db);
    int threads_per_block_val = 32; 
    dim3 threads_per_block_dim(threads_per_block_val); 
    dim3 block(32, 1, 1);
    dim3 grid(num_sequences_in_db, 1, 1);            // one block per DB sequence

    // round query length up to the nearest multiple of 32 for template dispatch
    int rounded_len = ((query_length + 31) >> 5) << 5;

    switch (rounded_len) {
        LAUNCH_CASE(32)   LAUNCH_CASE(64)   LAUNCH_CASE(96)
        LAUNCH_CASE(128)  LAUNCH_CASE(160)  LAUNCH_CASE(192)
        LAUNCH_CASE(224)  LAUNCH_CASE(256)  LAUNCH_CASE(288)
        LAUNCH_CASE(320)  LAUNCH_CASE(352)  LAUNCH_CASE(384)
        LAUNCH_CASE(416)  LAUNCH_CASE(448)  LAUNCH_CASE(480)
        LAUNCH_CASE(512)  LAUNCH_CASE(544)  LAUNCH_CASE(576)
        LAUNCH_CASE(608)  LAUNCH_CASE(640)  LAUNCH_CASE(672)
        LAUNCH_CASE(704)  LAUNCH_CASE(736)  LAUNCH_CASE(768)
        LAUNCH_CASE(800)  LAUNCH_CASE(832)  LAUNCH_CASE(864)
        LAUNCH_CASE(896)  LAUNCH_CASE(928)  LAUNCH_CASE(960)
        LAUNCH_CASE(992)  LAUNCH_CASE(1024)
    }
}

extern "C"
void launch_sw_cuda_affine(
        const uint8_t* query_seq_indices_ptr,
        int            query_length,
        int*           good_idx,
        const uint8_t* ascii,
        int*           lengths,
        int            num_sequences_in_db, // Actual number of DB sequences
        int*           output_scores_ptr,
        int            gap_open,
        int            gap_extend,
        cudaStream_t   stream)
{
    dim3 num_blocks(num_sequences_in_db);
    int threads_per_block_val = 32;
    dim3 threads_per_block_dim(threads_per_block_val);
    dim3 block(32, 1, 1);
    dim3 grid(num_sequences_in_db, 1, 1);
    int rounded_len = ((query_length + 31) >> 5) << 5;
    switch (rounded_len) {
        case 32:   sw_kernel_affine<32, 32>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 64:   sw_kernel_affine<32, 64>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 96:   sw_kernel_affine<32, 96>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 128:  sw_kernel_affine<32, 128> <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 160:  sw_kernel_affine<32, 160> <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 192:  sw_kernel_affine<32, 192> <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 224:  sw_kernel_affine<32, 224> <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 256:  sw_kernel_affine<32, 256> <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 288:  sw_kernel_affine<32, 288>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 320:  sw_kernel_affine<32, 320>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 352:  sw_kernel_affine<32, 352>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 384:  sw_kernel_affine<32, 384>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 416:  sw_kernel_affine<32, 416>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 448:  sw_kernel_affine<32, 448>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 480:  sw_kernel_affine<32, 480>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 512:  sw_kernel_affine<32, 512>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 544:  sw_kernel_affine<32, 544>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 576:  sw_kernel_affine<32, 576>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 608:  sw_kernel_affine<32, 608>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 640:  sw_kernel_affine<32, 640>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 672:  sw_kernel_affine<32, 672>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 704:  sw_kernel_affine<32, 704>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 736:  sw_kernel_affine<32, 736>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 768:  sw_kernel_affine<32, 768>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 800:  sw_kernel_affine<32, 800>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 832:  sw_kernel_affine<32, 832>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 864:  sw_kernel_affine<32, 864>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 896:  sw_kernel_affine<32, 896>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 928:  sw_kernel_affine<32, 928>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 960:  sw_kernel_affine<32, 960>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 992:  sw_kernel_affine<32, 992>  <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
        case 1024: sw_kernel_affine<32, 1024> <<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); break;
    }
}
