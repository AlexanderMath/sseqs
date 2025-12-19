#include <cuda.h>
#include <cuda_runtime.h>
#include <stdint.h> // For uint8_t
#include <cstdio>        
#include <cuda_fp16.h>

// Debug flag: set to 1 to enable debug printing in kernels
#define DEBUG 0

// this makes kernel 30% slower
#define PRINT_DEBUG 0 

// ============================================================================
// COMPILE FLAGS - Enable/disable individual kernels for faster compilation
// ============================================================================
#define COMPILE_LINEAR         0  // sw_kernel_linear (FP16) (screening)
#define COMPILE_AFFINE         1  // sw_kernel_affine (FP16)
#define COMPILE_AFFINE_INT8    0  // sw_kernel_affine_int8 (screening)
#define COMPILE_BACKTRACK      0  // sw_affine_backtrack_kernel
#define COMPILE_PROFILE        1  // sw_profile_kernel

// DEV MODE: Only compile 4 kernel sizes for fast iteration (5s compile vs 30s)
#define DEV_MODE               0  // Set to 0 for production (all sizes)
// ============================================================================

#define NUM_AMINO_ACIDS_CUDA 21  // 20 amino acids + 1 padding character
#define NUM_AMINO_ACIDS_CUDA_2 441  // 21*21

// Combined 3D BLOSUM packed in 4-bit: [qchar][db0][db1] -> uint8_t with both scores
// Low nibble = blosum[qchar][db0], high nibble = blosum[qchar][db1]
// 20*20*20 = 8000 bytes
__device__ uint8_t blosum62_3d_packed[NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA];

// Combined BLOSUM: [qchar][db0][db1] -> half2(blosum[q][d0], blosum[q][d1])
// 20*20*20 = 8000 half2 values = 32KB - fits in constant memory!
__device__ __constant__ __half2 blosum62_combined_cuda_constant[NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA];
__device__ __half2 blosum62_combined_cuda_global[NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA];

// Int8 BLOSUM for screening: [qchar][db0*20 + db1] -> uint16 packing 2x int8 (s0, s1)
// Layout: (int8_t s0, int8_t s1) packed as uint16 = (s1 << 8) | (s0 & 0xFF)
// 20 * 20 * 20 = 8000 uint16 values = 16KB
__device__ uint16_t blosum62_int8_combined_global[NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA];

// Simple 2D BLOSUM as half (not int8_t) - 800 bytes
__device__ half blosum62_as_half_global[NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA];

// 4-bit packed BLOSUM62 UPPER TRIANGLE (symmetric, 210 values -> 105 bytes, fits in 27 registers as uint32_t)
// Values offset by +4 to fit in 4 bits (0-15), subtract 4 after unpacking
// Stored as upper triangle: blosum[i][j] where i <= j
// Layout: [0,0], [0,1], [0,2], ..., [0,19], [1,1], [1,2], ..., [1,19], [2,2], ..., [19,19]
__device__ __constant__ uint32_t blosum62_upper_triangle_packed[27] = {
    // 210 values, 27 uint32_t - computed from BLOSUM62 matrix
    0x43342238, 0x53233332, 0x12494214, 0x36214245, 0xa1213321, 
    0x11544415, 0x20452124, 0x133641a1, 0x03431130, 0x31101d11, 
    0x23312313, 0x21426932, 0x32343145, 0x25114292, 0xa2213431, 
    0x42112002, 0x311c1122, 0x16223232, 0x32145168, 0x14628731, 
    0x13953232, 0x49221343, 0x0a533332, 0x33b37522, 0x22158210, 
    0xb16f4229, 0x00000083
};

// Row offsets for upper triangle indexing (precomputed)
__device__ __constant__ int blosum_row_offsets[20] = {
    0, 20, 39, 57, 74, 90, 105, 119, 132, 144, 
    155, 165, 174, 182, 189, 195, 200, 204, 207, 209
};

// Device function to get upper triangle index (OPTIMIZED)
__device__ __forceinline__ int blosum_upper_idx(int i, int j) {
    if (i > j) { int tmp = i; i = j; j = tmp; }
    return blosum_row_offsets[i] + (j - i);
}

// Device function to lookup from packed upper triangle (OPTIMIZED)
__device__ __forceinline__ int8_t blosum_lookup_packed(const uint32_t* packed, int qchar, int dbchar) {
    int idx = blosum_upper_idx(qchar, dbchar);
    int reg_idx = idx >> 3;  // Divide by 8 (faster than /)
    int bit_idx = (idx & 7) << 2;  // Modulo 8, then * 4 (faster than % and *)
    uint32_t val = (packed[reg_idx] >> bit_idx) & 0xF;
    return (int8_t)(val) - 4;  // Offset back to signed (-4 to +11)
}

// 4-bit packed BLOSUM62 lower triangle (symmetric, 210 values -> 105 bytes)
// Values offset by +4 to fit in 4 bits (0-15), subtract 4 after unpacking
__device__ __constant__ unsigned char blosum62_packed[105] = {
    0x34, 0x4e, 0xe2, 0x66, 0xe2, 0x74, 0xc0, 0xe0, 0xe0, 0x98, 0x05, 0x91, 0x08, 0x62, 0x63, 0x02,
    0x56, 0x03, 0x01, 0x06, 0xce, 0x66, 0x7d, 0x2c, 0xe1, 0x0c, 0xe3, 0xe1, 0x14, 0x09, 0xc9, 0x66,
    0xce, 0x21, 0x10, 0xce, 0xce, 0x61, 0x0e, 0x10, 0xce, 0xec, 0x15, 0xe6, 0x16, 0x84, 0xc0, 0xc2,
    0x41, 0xe5, 0x00, 0xc0, 0xce, 0xc0, 0xc1, 0xe1, 0xce, 0x01, 0xe5, 0x00, 0x00, 0xce, 0x41, 0xca,
    0xe0, 0xe5, 0x01, 0x81, 0xc0, 0xc1, 0x0c, 0xe3, 0xe1, 0xe1, 0xc1, 0xc1, 0xe2, 0xe0, 0x04, 0x05,
    0x0c, 0x0c, 0x04, 0xe5, 0x00, 0x0c, 0x0c, 0x06, 0xc0, 0x0b, 0xe0, 0x00, 0x34, 0x0c, 0x0c, 0xe5,
    0xe1, 0xe1, 0xf2, 0x65, 0x0c, 0xe2, 0xe0, 0x26
};

__device__ const int8_t blosum62_matrix_cuda_global[NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA] = {
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

// 4-bit packed BLOSUM62 (200 bytes): each byte stores 2 consecutive values
// Low nibble = even index, high nibble = odd index, add 4 before packing
__device__ const uint8_t blosum62_4bit_packed[200] = {
    0x38, 0x22, 0x34, 0x43, 0x32, 0x33, 0x23, 0x53, 0x14, 0x42, 
    0x93, 0x24, 0x51, 0x24, 0x14, 0x62, 0x13, 0x32, 0x13, 0x12, 
    0x42, 0x5a, 0x41, 0x44, 0x15, 0x41, 0x12, 0x52, 0x04, 0x12, 
    0x22, 0xa5, 0x41, 0x36, 0x13, 0x30, 0x11, 0x43, 0x03, 0x11, 
    0x14, 0x11, 0x1d, 0x10, 0x31, 0x13, 0x23, 0x31, 0x23, 0x32, 
    0x53, 0x44, 0x91, 0x26, 0x14, 0x52, 0x14, 0x43, 0x23, 0x23, 
    0x43, 0x64, 0x60, 0x29, 0x14, 0x51, 0x12, 0x43, 0x13, 0x22, 
    0x24, 0x34, 0x21, 0xa2, 0x02, 0x20, 0x11, 0x42, 0x22, 0x11, 
    0x42, 0x35, 0x41, 0x24, 0x1c, 0x31, 0x32, 0x32, 0x22, 0x16, 
    0x13, 0x11, 0x13, 0x01, 0x81, 0x16, 0x45, 0x21, 0x13, 0x73, 
    0x23, 0x01, 0x23, 0x01, 0x61, 0x28, 0x46, 0x21, 0x23, 0x53, 
    0x63, 0x34, 0x51, 0x25, 0x13, 0x92, 0x13, 0x43, 0x13, 0x22, 
    0x33, 0x12, 0x43, 0x12, 0x52, 0x36, 0x49, 0x32, 0x33, 0x53, 
    0x12, 0x11, 0x12, 0x11, 0x43, 0x14, 0xa4, 0x20, 0x52, 0x37, 
    0x23, 0x32, 0x31, 0x23, 0x12, 0x31, 0x02, 0x3b, 0x03, 0x21, 
    0x35, 0x45, 0x43, 0x44, 0x23, 0x42, 0x23, 0x83, 0x15, 0x22, 
    0x34, 0x34, 0x33, 0x23, 0x32, 0x33, 0x23, 0x53, 0x29, 0x42, 
    0x11, 0x00, 0x22, 0x21, 0x12, 0x12, 0x53, 0x10, 0xf2, 0x16, 
    0x22, 0x12, 0x32, 0x12, 0x36, 0x23, 0x73, 0x21, 0x62, 0x3b, 
    0x14, 0x11, 0x23, 0x12, 0x71, 0x25, 0x35, 0x22, 0x14, 0x83
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


// Include linear gap kernel
#include "sw_linear.cu"
#include "sw_affine.cu"
#include "sw_affine_uint8.cu"
#include "sw_profile.cu"

// Helper functions for half2 (Ada Lovelace optimized)

__device__ __forceinline__ int16_t half2int(half h) {
    return (int16_t)__half2short_rn(h);
}

__device__ __forceinline__ half int2half(int16_t i) {
    return __short2half_rn(i);
}

// Precomputed triangular numbers for unpacking: col*(col+1)/2 for col=0..19
__device__ __constant__ int tril_offsets[20] = { 0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66, 78, 91, 105, 120, 136, 153, 171, 190 };

// Simple unpacker for 4-bit packed BLOSUM from shared memory (symmetric matrix)
__device__ __forceinline__ int8_t unpack_blosum_shared(int row, int col, const unsigned char* packed) {
    // Ensure we access lower triangle (symmetric)
    if (row > col) {
        int tmp = row;
        row = col;
        col = tmp;
    }
    
    // Linear index in lower triangle using lookup table
    int linear_idx = tril_offsets[col] + row;
    int byte_idx = linear_idx >> 1;
    unsigned char packed_byte = packed[byte_idx];
    
    int8_t value;
    if (linear_idx & 1) {
        value = (int8_t)((packed_byte >> 4) & 0x0F) - 4;
    } else {
        value = (int8_t)(packed_byte & 0x0F) - 4;
    }
    
    return value;
}


// Unpack 4-bit BLOSUM value from packed array in registers
// Given (row, col) where row <= col (lower triangle), returns BLOSUM score
__device__ __forceinline__ int8_t unpack_blosum_4bit_from_regs(
    int row, int col,
    uint4 r0, uint4 r1, uint4 r2, uint4 r3, uint4 r4, uint4 r5, uint4 r6) {
    
    // Ensure we access lower triangle (symmetric matrix)
    if (row > col) {
        int tmp = row;
        row = col;
        col = tmp;
    }
    
    // Linear index in lower triangle
    int linear_idx = col * (col + 1) / 2 + row;
    
    // Two values per byte
    int byte_idx = linear_idx >> 1;
    
    // Extract byte from appropriate register
    unsigned char packed_byte;
    if (byte_idx < 16) {
        packed_byte = reinterpret_cast<const unsigned char*>(&r0)[byte_idx];
    } else if (byte_idx < 32) {
        packed_byte = reinterpret_cast<const unsigned char*>(&r1)[byte_idx - 16];
    } else if (byte_idx < 48) {
        packed_byte = reinterpret_cast<const unsigned char*>(&r2)[byte_idx - 32];
    } else if (byte_idx < 64) {
        packed_byte = reinterpret_cast<const unsigned char*>(&r3)[byte_idx - 48];
    } else if (byte_idx < 80) {
        packed_byte = reinterpret_cast<const unsigned char*>(&r4)[byte_idx - 64];
    } else if (byte_idx < 96) {
        packed_byte = reinterpret_cast<const unsigned char*>(&r5)[byte_idx - 80];
    } else {
        packed_byte = reinterpret_cast<const unsigned char*>(&r6)[byte_idx - 96];
    }
    
    // Extract 4-bit value
    int8_t value;
    if (linear_idx & 1) {
        value = (int8_t)((packed_byte >> 4) & 0x0F) - 4;
    } else {
        value = (int8_t)(packed_byte & 0x0F) - 4;
    }
    
    return value;
}

// ============================================================================
// Utility Kernels
// ============================================================================

// Simple kernel to zero out the output scores array
__global__ void zero_scores_kernel(int* out_scores, int num_sequences) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < num_sequences) {
        out_scores[idx] = 0;
    }
}


// Kernel to build combined BLOSUM [db0][db1][qchar] -> half2(blosum[q][d0], blosum[q][d1])
// NOTE: Dimension order changed to [db0][db1][qchar] to avoid bank conflicts!
// Also builds simple 2D BLOSUM as half
__global__ void build_combined_blosum_kernel() {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA;
    
    // Build simple 2D BLOSUM as half (first 400 threads)
    if (idx < NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA) {
        blosum62_as_half_global[idx] = int2half((int16_t)blosum62_matrix_cuda_global[idx]);
    }
    
    // Build combined 3D BLOSUM (half2 version - 32KB)
    // NEW LAYOUT: [db0][db1][qchar] instead of [qchar][db0][db1]
    if (idx < total) {
        int db0 = idx / (NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA);
        int remainder = idx % (NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA);
        int db1 = remainder / NUM_AMINO_ACIDS_CUDA;
        int qchar = remainder % NUM_AMINO_ACIDS_CUDA;
        
        int8_t score0, score1;
        
        // Padding character (index 20) gets large negative score to prevent accumulation
        if (db0 == 20 || qchar == 20) {
            score0 = -127;
        } else {
            score0 = blosum62_matrix_cuda_global[qchar * 20 + db0];  // Use 20, not NUM_AMINO_ACIDS_CUDA
        }
        
        if (db1 == 20 || qchar == 20) {
            score1 = -127;
        } else {
            score1 = blosum62_matrix_cuda_global[qchar * 20 + db1];
        }
        
        blosum62_combined_cuda_global[idx] = __halves2half2(int2half(score0), int2half(score1));
        
        // Build int8 combined version (16KB) - stores 2x int8_t scores as uint16
        uint16_t int8_packed = ((uint16_t)(score1 & 0xFF) << 8) | (uint16_t)(score0 & 0xFF);
        blosum62_int8_combined_global[idx] = int8_packed;
        
        // Build 4-bit packed version (8KB) - low nibble=score0, high nibble=score1
        uint8_t packed0 = (uint8_t)(score0 + 4);  // Offset to 0-15
        uint8_t packed1 = (uint8_t)(score1 + 4);
        blosum62_3d_packed[idx] = (packed1 << 4) | packed0;
    }
}

// Host function to initialize combined BLOSUM (call once at startup)
extern "C" void initialize_combined_blosum() {
    int total = NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA;
    int threads = 256;
    int blocks = (total + threads - 1) / threads;
    build_combined_blosum_kernel<<<blocks, threads>>>();
    cudaDeviceSynchronize();
    
    // Copy to constant memory for faster access
    cudaMemcpyToSymbol(blosum62_combined_cuda_constant, blosum62_combined_cuda_global, 
                       total * sizeof(__half2));
    cudaDeviceSynchronize();
}



#define LAUNCH_CASE(LEN)                                                         \
    case LEN:                                                                    \
        sw_kernel<32, LEN><<<grid, block, 0, stream>>>(query_seq_indices_ptr, query_length, good_idx, ascii, lengths, num_sequences_in_db, output_scores_ptr, gap_val);                                                        \
        break;



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
        int            db_seq_length,  // NEW: sequence length passed from Python
        cudaStream_t   stream)
{
    // Initialize 3D packed BLOSUM once per device
    static bool blosum_initialized[8] = {false}; // Support up to 8 GPUs
    int device;
    cudaGetDevice(&device);
    if (!blosum_initialized[device]) {
        int total = NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA;
        int threads = 256;
        int blocks = (total + threads - 1) / threads;
        build_combined_blosum_kernel<<<blocks, threads>>>();
        cudaDeviceSynchronize();
        blosum_initialized[device] = true;
    }
    
    // Zero out output scores to avoid garbage from unprocessed sequences
    {
        int threads = 256;
        int blocks = (num_sequences_in_db + threads - 1) / threads;
        zero_scores_kernel<<<blocks, threads, 0, stream>>>(output_scores_ptr, num_sequences_in_db);
    }
    
    // 256 threads per block (8 warps), processing 64 sequences per block
    constexpr int THREADS_PER_BLOCK = 256;
    
    // Use the sequence length passed from Python
    // printf("DEBUG: db_seq_length=%d (from Python)\n", db_seq_length);
    
    // Round up to nearest multiple of 8, then find closest supported length
    int rounded = ((db_seq_length + 7) / 8) * 8;

    // Round up to nearest multiple of 32
    int dispatch_length = ((db_seq_length + 31) / 32) * 32;

    // Clamp to valid range (kernel will handle individual sequences via length0/length1)
    //if (dispatch_length < 32) dispatch_length = 32;
    //if (dispatch_length > 1024) dispatch_length = 1024;

    // Macro to launch kernel with given SUB_WARP_SIZE and TPT
    #define LAUNCH_KERNEL(SUB_WARP, TPT_VAL) \
        do { \
            constexpr int SUB_WARP_SIZE = SUB_WARP; \
            constexpr int TPT = TPT_VAL; \
            constexpr int SEQS_PER_BLOCK = (THREADS_PER_BLOCK / 32) * (32 / SUB_WARP_SIZE) * 2; \
            dim3 block(THREADS_PER_BLOCK, 1, 1); \
            dim3 grid((num_sequences_in_db + SEQS_PER_BLOCK - 1) / SEQS_PER_BLOCK, 1, 1); \
            sw_kernel_affine<THREADS_PER_BLOCK, SUB_WARP_SIZE, TPT> \
                <<<grid, block, 0, stream>>>( \
                    query_seq_indices_ptr, query_length, good_idx, ascii, lengths, \
                    num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); \
        } while(0)
    
    // Dispatch based on rounded-up database sequence length (SUB_WARP_SIZE * TPT >= db_seq_length)
    switch (dispatch_length) {
        case 32: LAUNCH_KERNEL(8, 4); break;  // ratio: 2.0
        case 64: LAUNCH_KERNEL(8, 8); break;  // ratio: 1.0
        case 96: LAUNCH_KERNEL(8, 12); break;  // ratio: 1.5
        case 128: LAUNCH_KERNEL(8, 16); break;  // ratio: 2.0
        case 160: LAUNCH_KERNEL(16, 10); break;  // ratio: 1.6
        case 192: LAUNCH_KERNEL(16, 12); break;  // ratio: 1.3
        case 224: LAUNCH_KERNEL(16, 14); break;  // ratio: 1.1
        case 256: LAUNCH_KERNEL(16, 16); break;  // ratio: 1.0
        case 288: LAUNCH_KERNEL(16, 18); break;  // ratio: 1.1
        case 320: LAUNCH_KERNEL(16, 20); break;  // ratio: 1.2
        case 352: LAUNCH_KERNEL(16, 22); break;  // ratio: 1.4
        case 384: LAUNCH_KERNEL(16, 24); break;  // ratio: 1.5
        case 416: LAUNCH_KERNEL(16, 26); break;  // ratio: 1.6
        case 448: LAUNCH_KERNEL(16, 28); break;  // ratio: 1.8
        case 480: LAUNCH_KERNEL(16, 30); break;  // ratio: 1.9
        case 512: LAUNCH_KERNEL(16, 32); break;  // ratio: 2.0
        case 544: LAUNCH_KERNEL(32, 17); break;  // ratio: 1.9
        case 576: LAUNCH_KERNEL(32, 18); break;  // ratio: 1.8
        case 608: LAUNCH_KERNEL(32, 19); break;  // ratio: 1.7
        case 640: LAUNCH_KERNEL(32, 20); break;  // ratio: 1.6
        case 672: LAUNCH_KERNEL(32, 21); break;  // ratio: 1.5
        case 704: LAUNCH_KERNEL(32, 22); break;  // ratio: 1.5
        case 736: LAUNCH_KERNEL(32, 23); break;  // ratio: 1.4
        case 768: LAUNCH_KERNEL(32, 24); break;  // ratio: 1.3
        case 800: LAUNCH_KERNEL(32, 25); break;  // ratio: 1.3
        case 832: LAUNCH_KERNEL(32, 26); break;  // ratio: 1.2
        case 864: LAUNCH_KERNEL(32, 27); break;  // ratio: 1.2
        case 896: LAUNCH_KERNEL(32, 28); break;  // ratio: 1.1
        case 928: LAUNCH_KERNEL(32, 29); break;  // ratio: 1.1
        case 960: LAUNCH_KERNEL(32, 30); break;  // ratio: 1.1
        case 992: LAUNCH_KERNEL(32, 31); break;  // ratio: 1.0
        case 1024: LAUNCH_KERNEL(32, 32); break;  // ratio: 1.0
        default:
            printf("Error: Unsupported dispatch_length %d\n", dispatch_length);
            return;
    }

        
    #undef LAUNCH_KERNEL
}



extern "C"
void launch_sw_cuda_profile(
        const int8_t*  pssm_ptr,        // PSSM matrix: (query_length x 20)
        int            query_length,
        int*           good_idx,
        const uint8_t* ascii,
        int*           lengths,
        int            num_sequences_in_db,
        int*           output_scores_ptr,
        int            gap_open,
        int            gap_extend,
        int            db_seq_length,
        cudaStream_t   stream)
{
    // No BLOSUM initialization needed - we use PSSM
    
    // Zero out output scores to avoid garbage from unprocessed sequences
    {
        int threads = 256;
        int blocks = (num_sequences_in_db + threads - 1) / threads;
        zero_scores_kernel<<<blocks, threads, 0, stream>>>(output_scores_ptr, num_sequences_in_db);
    }
    
    // 256 threads per block (8 warps), processing 64 sequences per block
    constexpr int THREADS_PER_BLOCK = 256;
    
    // Round up to nearest multiple of 32
    int dispatch_length = ((db_seq_length + 31) / 32) * 32;

    // Macro to launch kernel with given SUB_WARP_SIZE and TPT
    #define LAUNCH_KERNEL(SUB_WARP, TPT_VAL) \
        do { \
            constexpr int SUB_WARP_SIZE = SUB_WARP; \
            constexpr int TPT = TPT_VAL; \
            constexpr int SEQS_PER_BLOCK = (THREADS_PER_BLOCK / 32) * (32 / SUB_WARP_SIZE) * 2; \
            dim3 block(THREADS_PER_BLOCK, 1, 1); \
            dim3 grid((num_sequences_in_db + SEQS_PER_BLOCK - 1) / SEQS_PER_BLOCK, 1, 1); \
            sw_kernel_profile<THREADS_PER_BLOCK, SUB_WARP_SIZE, TPT> \
                <<<grid, block, 0, stream>>>( \
                    pssm_ptr, query_length, good_idx, ascii, lengths, \
                    num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); \
        } while(0)
    
    // Dispatch based on rounded-up database sequence length (SUB_WARP_SIZE * TPT >= db_seq_length)
    switch (dispatch_length) {
        case 32: LAUNCH_KERNEL(8, 4); break;  // ratio: 2.0
        case 64: LAUNCH_KERNEL(8, 8); break;  // ratio: 1.0
        case 96: LAUNCH_KERNEL(8, 12); break;  // ratio: 1.5
        case 128: LAUNCH_KERNEL(8, 16); break;  // ratio: 2.0
        case 160: LAUNCH_KERNEL(16, 10); break;  // ratio: 1.6
        case 192: LAUNCH_KERNEL(16, 12); break;  // ratio: 1.3
        case 224: LAUNCH_KERNEL(16, 14); break;  // ratio: 1.1
        case 256: LAUNCH_KERNEL(16, 16); break;  // ratio: 1.0
        case 288: LAUNCH_KERNEL(16, 18); break;  // ratio: 1.1
        case 320: LAUNCH_KERNEL(16, 20); break;  // ratio: 1.2
        case 352: LAUNCH_KERNEL(16, 22); break;  // ratio: 1.4
        case 384: LAUNCH_KERNEL(16, 24); break;  // ratio: 1.5
        case 416: LAUNCH_KERNEL(16, 26); break;  // ratio: 1.6
        case 448: LAUNCH_KERNEL(16, 28); break;  // ratio: 1.8
        case 480: LAUNCH_KERNEL(16, 30); break;  // ratio: 1.9
        case 512: LAUNCH_KERNEL(16, 32); break;  // ratio: 2.0
        case 544: LAUNCH_KERNEL(32, 17); break;  // ratio: 1.9
        case 576: LAUNCH_KERNEL(32, 18); break;  // ratio: 1.8
        case 608: LAUNCH_KERNEL(32, 19); break;  // ratio: 1.7
        case 640: LAUNCH_KERNEL(32, 20); break;  // ratio: 1.6
        case 672: LAUNCH_KERNEL(32, 21); break;  // ratio: 1.5
        case 704: LAUNCH_KERNEL(32, 22); break;  // ratio: 1.5
        case 736: LAUNCH_KERNEL(32, 23); break;  // ratio: 1.4
        case 768: LAUNCH_KERNEL(32, 24); break;  // ratio: 1.3
        case 800: LAUNCH_KERNEL(32, 25); break;  // ratio: 1.3
        case 832: LAUNCH_KERNEL(32, 26); break;  // ratio: 1.2
        case 864: LAUNCH_KERNEL(32, 27); break;  // ratio: 1.2
        case 896: LAUNCH_KERNEL(32, 28); break;  // ratio: 1.1
        case 928: LAUNCH_KERNEL(32, 29); break;  // ratio: 1.1
        case 960: LAUNCH_KERNEL(32, 30); break;  // ratio: 1.1
        case 992: LAUNCH_KERNEL(32, 31); break;  // ratio: 1.0
        case 1024: LAUNCH_KERNEL(32, 32); break;  // ratio: 1.0
        default:
            printf("Error: Unsupported dispatch_length %d\n", dispatch_length);
            return;
    }

        
    #undef LAUNCH_KERNEL
}

// ============================================================================
// Linear Gap Penalty Launcher (for screening/prefilter)
// ============================================================================
extern "C"
void launch_sw_cuda_linear(
        const uint8_t* query_seq_indices_ptr,
        int            query_length,
        int*           good_idx,
        const uint8_t* ascii,
        int*           lengths,
        int            num_sequences_in_db,
        int*           output_scores_ptr,
        int            gap_penalty,       // Single gap penalty (linear)
        int            db_seq_length,     // Sequence length passed from Python
        cudaStream_t   stream)
{
#if COMPILE_LINEAR
    // Initialize 3D packed BLOSUM once
    static bool blosum_initialized = false;
    if (!blosum_initialized) {
        int total = NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA;
        int threads = 256;
        int blocks = (total + threads - 1) / threads;
        build_combined_blosum_kernel<<<blocks, threads>>>();
        cudaDeviceSynchronize();
        
        // Copy to constant memory for linear kernel
        cudaMemcpyToSymbol(blosum62_combined_cuda_constant, blosum62_combined_cuda_global, 
                           total * sizeof(__half2));
        cudaDeviceSynchronize();
        
        blosum_initialized = true;
    }
    
    // Zero out output scores to avoid garbage from unprocessed sequences
    {
        int threads = 256;
        int blocks = (num_sequences_in_db + threads - 1) / threads;
        zero_scores_kernel<<<blocks, threads, 0, stream>>>(output_scores_ptr, num_sequences_in_db);
    }
    
    // 256 threads per block (8 warps), processing sequences per block
    constexpr int THREADS_PER_BLOCK = 256;
    
    // Round up to nearest multiple of 32
    int dispatch_length = ((db_seq_length + 31) / 32) * 32;

    // Macro to launch kernel with given SUB_WARP_SIZE and TPT
    #define LAUNCH_KERNEL(SUB_WARP, TPT_VAL) \
        do { \
            constexpr int SUB_WARP_SIZE = SUB_WARP; \
            constexpr int TPT = TPT_VAL; \
            constexpr int SEQS_PER_BLOCK = (THREADS_PER_BLOCK / 32) * (32 / SUB_WARP_SIZE) * 2; \
            dim3 block(THREADS_PER_BLOCK, 1, 1); \
            dim3 grid((num_sequences_in_db + SEQS_PER_BLOCK - 1) / SEQS_PER_BLOCK, 1, 1); \
            sw_kernel_linear<THREADS_PER_BLOCK, SUB_WARP_SIZE, TPT> \
                <<<grid, block, 0, stream>>>( \
                    query_seq_indices_ptr, query_length, good_idx, ascii, lengths, \
                    num_sequences_in_db, output_scores_ptr, gap_penalty); \
        } while(0)
    
    // Dispatch based on rounded-up database sequence length (SUB_WARP_SIZE * TPT >= db_seq_length)
    // LINEAR GAP: Using 2x higher TPT than affine due to ~50% lower register pressure
    
    #if DEV_MODE
    // DEV MODE: Only compile 4 sizes for fast iteration
    switch (dispatch_length) {
        case 32:
        case 64:
        case 96:
        case 128: LAUNCH_KERNEL(8, 16); break;
        case 160:
        case 192:
        case 224:
        case 256: LAUNCH_KERNEL(8, 32); break;
        case 288:
        case 320:
        case 352:
        case 384:
        case 416:
        case 448:
        case 480:
        case 512: LAUNCH_KERNEL(16, 32); break;
        default:
            // 544-1024 and beyond
            LAUNCH_KERNEL(16, 64);
            break;
    }
    #else
    // PRODUCTION MODE: Compile optimized size for each length
    switch (dispatch_length) {
        case 32: LAUNCH_KERNEL(8, 4); break;
        case 64: LAUNCH_KERNEL(8, 8); break;
        case 96: LAUNCH_KERNEL(8, 12); break;
        case 128: LAUNCH_KERNEL(8, 16); break;
        case 160: LAUNCH_KERNEL(8, 20); break;
        case 192: LAUNCH_KERNEL(8, 24); break;
        case 224: LAUNCH_KERNEL(8, 28); break;
        case 256: LAUNCH_KERNEL(8, 32); break;
        case 288: LAUNCH_KERNEL(16, 18); break;
        case 320: LAUNCH_KERNEL(16, 20); break;
        case 352: LAUNCH_KERNEL(16, 22); break;
        case 384: LAUNCH_KERNEL(16, 24); break;
        case 416: LAUNCH_KERNEL(16, 26); break;
        case 448: LAUNCH_KERNEL(16, 28); break;
        case 480: LAUNCH_KERNEL(16, 30); break;
        case 512: LAUNCH_KERNEL(16, 32); break;
        case 544: LAUNCH_KERNEL(16, 34); break;
        case 576: LAUNCH_KERNEL(16, 36); break;
        case 608: LAUNCH_KERNEL(16, 38); break;
        case 640: LAUNCH_KERNEL(16, 40); break;
        case 672: LAUNCH_KERNEL(16, 42); break;
        case 704: LAUNCH_KERNEL(16, 44); break;
        case 736: LAUNCH_KERNEL(16, 46); break;
        case 768: LAUNCH_KERNEL(16, 48); break;
        case 800: LAUNCH_KERNEL(16, 50); break;
        case 832: LAUNCH_KERNEL(16, 52); break;
        case 864: LAUNCH_KERNEL(16, 54); break;
        case 896: LAUNCH_KERNEL(16, 56); break;
        case 928: LAUNCH_KERNEL(16, 58); break;
        case 960: LAUNCH_KERNEL(16, 60); break;
        case 992: LAUNCH_KERNEL(16, 62); break;
        case 1024: LAUNCH_KERNEL(16, 64); break;
        default:
            printf("Error: Unsupported dispatch_length %d\n", dispatch_length);
            return;
    }
    #endif  // DEV_MODE
        
    #undef LAUNCH_KERNEL
#else
    printf("ERROR: COMPILE_LINEAR is disabled!\n");
#endif
}



extern "C"
void launch_sw_cuda_affine_uint8(
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
    // Initialize 3D packed BLOSUM once
    static bool blosum_initialized = false;
    if (!blosum_initialized) {
        int total = NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA;
        int threads = 256;
        int blocks = (total + threads - 1) / threads;
        build_combined_blosum_kernel<<<blocks, threads>>>();
        cudaDeviceSynchronize();
        blosum_initialized = true;
    }
    
    // Zero out output scores to avoid garbage from unprocessed sequences
    {
        int threads = 256;
        int blocks = (num_sequences_in_db + threads - 1) / threads;
        zero_scores_kernel<<<blocks, threads, 0, stream>>>(output_scores_ptr, num_sequences_in_db);
    }
    
    // 256 threads per block (8 warps), processing 64 sequences per block
    constexpr int THREADS_PER_BLOCK = 256;
    
    // Compute database sequence length from first sequence
    int h_starts[2];
    cudaMemcpy(h_starts, lengths, 2 * sizeof(int), cudaMemcpyDeviceToHost);
    int db_seq_length = h_starts[1] - h_starts[0] - 1;
    
    // Macro to launch kernel with given SUB_WARP_SIZE and TPT
    #define LAUNCH_KERNEL(SUB_WARP, TPT_VAL) \
        do { \
            constexpr int SUB_WARP_SIZE = SUB_WARP; \
            constexpr int TPT = TPT_VAL; \
            constexpr int SEQS_PER_BLOCK = (THREADS_PER_BLOCK / 32) * (32 / SUB_WARP_SIZE) * 4; \
            dim3 block(THREADS_PER_BLOCK, 1, 1); \
            dim3 grid((num_sequences_in_db + SEQS_PER_BLOCK - 1) / SEQS_PER_BLOCK, 1, 1); \
            sw_kernel_affine_int8<THREADS_PER_BLOCK, SUB_WARP_SIZE, TPT> \
                <<<grid, block, 0, stream>>>( \
                    query_seq_indices_ptr, query_length, good_idx, ascii, lengths, \
                    num_sequences_in_db, output_scores_ptr, gap_open, gap_extend); \
        } while(0)
    
    // Dispatch based on database sequence length (SUB_WARP_SIZE * TPT = db_seq_length)
    switch (db_seq_length) {
        case 128: LAUNCH_KERNEL(8, 16);  break;
        case 256: LAUNCH_KERNEL(16, 16); break;
        case 512: LAUNCH_KERNEL(32, 16); break;
        case 1024: LAUNCH_KERNEL(32, 32); break;
        default:
            printf("Error: Unsupported db_seq_length %d\n", db_seq_length);
            return;
    }
    
    #undef LAUNCH_KERNEL
}

extern "C"
void launch_sw_cuda_affine_antidiag(
        const uint8_t* query_seq_indices_ptr,
        int            query_length,
        int*           good_idx,
        const uint8_t* ascii,
        int*           lengths,
        int            num_sequences_in_db,
        int*           output_scores_ptr,
        int            gap_open,
        int            gap_extend,
        cudaStream_t   stream) {

            return;
}

// ================================================================================================
// BACKTRACKING VERSION - 
// ================================================================================================
#include "sw_backtrack_vram.cu"

extern "C"
void launch_sw_cuda_affine_backtrack(
        const uint8_t* query_seq_indices_ptr,
        int            query_length,
        int*           good_idx,
        const uint8_t* ascii,
        int*           lengths,
        int            num_sequences_in_db,
        int*           output_scores_ptr,
        int*           output_alignment_lens_ptr,
        int*           output_end_i_ptr,
        int*           output_end_j_ptr,
        char*          output_alignment_ops_ptr,
        int            max_align_len,
        int            gap_open,
        int            gap_extend,
        cudaStream_t   stream)
{
    // Initialize 3D packed BLOSUM once per device
    static bool blosum_initialized[8] = {false}; // Support up to 8 GPUs
    int device;
    cudaGetDevice(&device);
    if (!blosum_initialized[device]) {
        int total = NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA * NUM_AMINO_ACIDS_CUDA;
        int threads = 256;
        int blocks = (total + threads - 1) / threads;
        build_combined_blosum_kernel<<<blocks, threads>>>();
        cudaDeviceSynchronize();
        blosum_initialized[device] = true;
    }
    
    // Zero out output scores to avoid garbage from unprocessed sequences
    {
        int threads = 256;
        int blocks = (num_sequences_in_db + threads - 1) / threads;
        zero_scores_kernel<<<blocks, threads, 0, stream>>>(output_scores_ptr, num_sequences_in_db);
    }
    
    // Allocate matrix buffers in global memory (instead of checkpoints)
    // Each sequence needs: [max_target_len][max_query_len] matrices for H, E, F
    // Using half2 (4 bytes each), storing packed values
    
    // 256 threads per block (8 warps), processing sequences per block
    constexpr int THREADS_PER_BLOCK = 256;
    
    // Compute actual database sequence length from first sequence and find max
    // Use heap allocation for large arrays (avoid stack overflow)
    int* h_starts = new int[num_sequences_in_db + 1];
    cudaMemcpy(h_starts, lengths, (num_sequences_in_db + 1) * sizeof(int), cudaMemcpyDeviceToHost);
    
    int max_target_len = 0;
    for (int i = 0; i < num_sequences_in_db; i++) {
        int seq_len = h_starts[i + 1] - h_starts[i] - 1;
        if (seq_len > max_target_len) max_target_len = seq_len;
    }
    
    int db_seq_length = h_starts[1] - h_starts[0] - 1;  // First sequence length for dispatch
    delete[] h_starts;  // Free heap memory
    
    int max_query_len = query_length;
    
    // Round up max_target_len to nearest 32 for better memory alignment
    max_target_len = ((max_target_len + 31) / 32) * 32;
    
    // Storage: direction pointers (1 byte per cell, contains 2x2-bit directions)
    // Matrix dimensions: [num_pairs][max_query_len][max_target_len]
    // Note: We use double-packing (2 seqs per pair), so allocate based on pairs not sequences
    int num_pairs = (num_sequences_in_db + 1) / 2;  // Round up for odd number of sequences
    size_t matrix_size_per_pair = (size_t)max_query_len * max_target_len;
    size_t matrix_bytes = num_pairs * matrix_size_per_pair * sizeof(uint8_t);  // 1 byte per cell!
    
    uint8_t* global_directions;
    cudaMalloc(&global_directions, matrix_bytes);
    
    // Round up to nearest multiple of 32
    int dispatch_length = ((db_seq_length + 31) / 32) * 32;

    // Macro to launch kernel with given SUB_WARP_SIZE and TPT
    #define LAUNCH_KERNEL(SUB_WARP, TPT_VAL) \
        do { \
            constexpr int SUB_WARP_SIZE = SUB_WARP; \
            constexpr int TPT = TPT_VAL; \
            constexpr int SEQS_PER_BLOCK = (THREADS_PER_BLOCK / 32) * (32 / SUB_WARP_SIZE) * 2; \
            dim3 block(THREADS_PER_BLOCK, 1, 1); \
            dim3 grid((num_sequences_in_db + SEQS_PER_BLOCK - 1) / SEQS_PER_BLOCK, 1, 1); \
            sw_kernel_affine_backtrack<THREADS_PER_BLOCK, SUB_WARP_SIZE, TPT> \
                <<<grid, block, 0, stream>>>( \
                    query_seq_indices_ptr, query_length, good_idx, ascii, lengths, \
                    num_sequences_in_db, output_scores_ptr, output_alignment_lens_ptr, \
                    output_end_i_ptr, output_end_j_ptr, output_alignment_ops_ptr, max_align_len, \
                    gap_open, gap_extend, global_directions, max_target_len, max_query_len); \
        } while(0)
    
    // Dispatch based on rounded-up database sequence length (SUB_WARP_SIZE * TPT >= db_seq_length)
    // DEV MODE: Only compile 4 sizes for fast iteration (~5s vs 88s compile time)
    #if DEV_MODE
    switch (dispatch_length) {
        case 32:
        case 64:
        case 96:
        case 128: LAUNCH_KERNEL(8, 16); break;
        case 160:
        case 192:
        case 224:
        case 256: LAUNCH_KERNEL(16, 16); break;
        case 288:
        case 320:
        case 352:
        case 384:
        case 416:
        case 448:
        case 480:
        case 512: LAUNCH_KERNEL(16, 32); break;
        default:
            // 544-1024 and beyond
            LAUNCH_KERNEL(32, 32);
            break;
    }
    #else
    // PRODUCTION MODE: Compile optimized kernel for each length
    switch (dispatch_length) {
        case 32: LAUNCH_KERNEL(8, 4); break;
        case 64: LAUNCH_KERNEL(8, 8); break;
        case 96: LAUNCH_KERNEL(8, 12); break;
        case 128: LAUNCH_KERNEL(8, 16); break;
        case 160: LAUNCH_KERNEL(16, 10); break;
        case 192: LAUNCH_KERNEL(16, 12); break;
        case 224: LAUNCH_KERNEL(16, 14); break;
        case 256: LAUNCH_KERNEL(16, 16); break;
        case 288: LAUNCH_KERNEL(16, 18); break;
        case 320: LAUNCH_KERNEL(16, 20); break;
        case 352: LAUNCH_KERNEL(16, 22); break;
        case 384: LAUNCH_KERNEL(16, 24); break;
        case 416: LAUNCH_KERNEL(16, 26); break;
        case 448: LAUNCH_KERNEL(16, 28); break;
        case 480: LAUNCH_KERNEL(16, 30); break;
        case 512: LAUNCH_KERNEL(16, 32); break;
        case 544: LAUNCH_KERNEL(32, 17); break;
        case 576: LAUNCH_KERNEL(32, 18); break;
        case 608: LAUNCH_KERNEL(32, 19); break;
        case 640: LAUNCH_KERNEL(32, 20); break;
        case 672: LAUNCH_KERNEL(32, 21); break;
        case 704: LAUNCH_KERNEL(32, 22); break;
        case 736: LAUNCH_KERNEL(32, 23); break;
        case 768: LAUNCH_KERNEL(32, 24); break;
        case 800: LAUNCH_KERNEL(32, 25); break;
        case 832: LAUNCH_KERNEL(32, 26); break;
        case 864: LAUNCH_KERNEL(32, 27); break;
        case 896: LAUNCH_KERNEL(32, 28); break;
        case 928: LAUNCH_KERNEL(32, 29); break;
        case 960: LAUNCH_KERNEL(32, 30); break;
        case 992: LAUNCH_KERNEL(32, 31); break;
        case 1024: LAUNCH_KERNEL(32, 32); break;
        default:
            printf("Error: Unsupported dispatch_length %d\n", dispatch_length);
            return;
    }
    #endif  // DEV_MODE
        
    #undef LAUNCH_KERNEL
    
    // Launch traceback kernel to reconstruct alignments
    {
        int threads = 256;
        int blocks = (num_sequences_in_db + threads - 1) / threads;
        sw_traceback_kernel<<<blocks, threads, 0, stream>>>(
            global_directions,
            output_scores_ptr,
            output_alignment_lens_ptr,
            output_end_i_ptr,
            output_end_j_ptr,
            output_alignment_ops_ptr,
            max_align_len,
            num_sequences_in_db,
            max_target_len,
            max_query_len
        );
    }
    
    // Free direction buffer
    cudaFree(global_directions);
}