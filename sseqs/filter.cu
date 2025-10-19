#include <cuda_runtime.h>
#include <cuda_fp16.h>
#include <cstdint>
#include <cfloat>
#include <cstdio>
#include "fft.cu"

#define M_PI 3.14159265358979323846
#define things_per_thread 16


template<int LOG_H, bool inverse = false>
__device__ inline void fft_inplace_warp(__half2* s)
{
    constexpr int  H      = 1 << LOG_H;   // total complex points
    constexpr int  TPT    = H / 32;       // how many each lane owns
    constexpr unsigned ALL = 0xffffffffu;

    const int lane = threadIdx.x & 31;    // 0 … 31
    const int base = lane * TPT;          // first index this lane owns

    /* -------- 1. load my TPT values into registers ---------------- */
    __half2 v[TPT];
#pragma unroll
    for (int i = 0; i < TPT; ++i)
        v[i] = s[base + i];

    /* -------- 2. radix-2 stages ---------------------------------- */
#pragma unroll
    for (int stage = 0; stage < LOG_H; ++stage)
    {
        const int m    = 1 << (stage + 1);   // 2,4,8,…
        const int half = m >> 1;             // 1,2,4,…
        //const float k2π_over_m = -2.0f * (float)M_PI / m;
        const float k2π_over_m = (inverse ? +2.0f : -2.0f) * (float)M_PI / m;


        if (half < TPT)
        {
            /* partner is another element **inside the same lane** */
#pragma unroll
            for (int i = 0; i < TPT; ++i)
            {
                int g = base + i;
                if ((g & half) == 0)
                {
                    int j = g & (half - 1);
                    float θ = k2π_over_m * j;
                    __half2 w = h2(::cosf(θ), ::sinf(θ));

                    __half2 a = v[i];
                    __half2 b = v[i ^ half];
                    __half2 t = cmul(b, w);

                    v[i]         = add(a, t);
                    v[i ^ half]  = sub(a, t);
                }
            }
        }
        else
        {
            /* partner lives in ANOTHER lane → one shuffle          */
            const int lane_delta = half / TPT;
#pragma unroll
            for (int i = 0; i < TPT; ++i)
            {
                int g = base + i;
                bool lower = (g & half) == 0;

                __half2 a = v[i];
                __half2 b = __shfl_xor_sync(ALL, a, lane_delta);

                int   j = g & (half - 1);
                float θ = k2π_over_m * j;
                __half2 w = h2(::cosf(θ), ::sinf(θ));

                if (lower) {
                    v[i] = add(a, cmul(b, w));      // a + b·w
                } else {
                    v[i] = sub(b, cmul(a, w));      // b – a·w
                }
            }
        }
    }

    /* -------- 3. store results back to shared memory -------------- */
#pragma unroll
    for (int i = 0; i < TPT; ++i)
        s[base + i] = v[i];
}

template<int LOG_N>
__global__ void precompute_query_rfft_kernel(
    __half2*               __restrict__ query_rfft_coeffs,    // (k, H+1) output
    const uint8_t*         __restrict__ query_reversed_indices,
    const float*           __restrict__ weight_LUTS_A_flat,
    int  k,            // number of ranks
    int  vocab_size,   // vocab_size_padded
    int  query_len)    // actual length of query sequence
{
    constexpr int  N      = 1 << LOG_N;       // time-domain length (power-of-2)
    constexpr int  H      = N >> 1;           // complex length after pairing (N/2)
    constexpr int  LOG_H  = LOG_N - 1;

    int bid = blockIdx.x;        // rank id (0 to k-1)
    int tid = threadIdx.x;       // 0 … H-1

    if (bid >= k || tid >= H) return;

    // Shared memory for this rank's query processing
    extern __shared__ char _shmem[];
    __half2* s_query_packed_time = (__half2*)_shmem;        // H elements

    /* A.1. Gather & Encode Query tokens for rank bid */
    {
        int idx0 = 2 * tid;
        float q_w0 = (idx0 < query_len)
                     ? weight_LUTS_A_flat[bid * vocab_size + query_reversed_indices[idx0]]
                     : 0.f;
        float q_w1 = (idx0 + 1 < query_len)
                     ? weight_LUTS_A_flat[bid * vocab_size + query_reversed_indices[idx0 + 1]]
                     : 0.f;
        s_query_packed_time[tid] = h2(q_w0, q_w1); 
    }
    __syncthreads();

    /* A.2. Forward FFT on H points for Query */
    {
        int rev = __brev(tid) >> (32 - LOG_H);
        if (tid < rev) {
            __half2 tmp_val = s_query_packed_time[tid];
            s_query_packed_time[tid] = s_query_packed_time[rev];
            s_query_packed_time[rev] = tmp_val;
        }
        __syncthreads();
        //fft_inplace<LOG_H>(s_query_packed_time, tid);
        fft_inplace_warp<LOG_H, false>(s_query_packed_time);
    }

    /* A.3. Unpack G_query[k] to Query_rfft[k] and store to global memory */
    __half2* output_row = query_rfft_coeffs + bid * (H + 1);
    
    if (tid == 0) {
        __half2 G0_q = s_query_packed_time[0];
        output_row[0] = h2(real(G0_q) + imag(G0_q), 0.f); 
        output_row[H] = h2(real(G0_q) - imag(G0_q), 0.f); 
    } else { // tid = 1 .. H-1
        int k_idx = tid;
        __half2 Gk_q_val = s_query_packed_time[k_idx];
        __half2 GH_minus_k_q_conj_val = conjh(s_query_packed_time[H - k_idx]);
        float theta = -2.0f * (float)M_PI * k_idx / N;
        __half2 WNk = h2(::cosf(theta), ::sinf(theta));
        __half2 term_sum  = add(Gk_q_val, GH_minus_k_q_conj_val);
        __half2 term_diff = sub(Gk_q_val, GH_minus_k_q_conj_val);
        __half2 j_WNk_term_diff = cmul(h2(0.f, -1.f), cmul(WNk, term_diff));
        output_row[k_idx] = h2(0.5f * real(add(term_sum, j_WNk_term_diff)), 
                              0.5f * imag(add(term_sum, j_WNk_term_diff)));
    }
}

extern "C"
void launch_warp_corr(
    __half* corr_ptr,
    const uint8_t* query_reversed_indices_ptr, // NEW
    const float*   weight_LUTS_A_ptr,          // NEW
    const float*   W_db_ptr,                   // for DB seqs (_weight_LUTS_B)
    const uint8_t* flat_db_ptr,
    const int32_t* offsets_db_ptr,
    const int32_t* lengths_db_ptr,
    float*         out_stats_ptr, // unused
    int B, int vocab_size, int k, int query_len, 
    int maxlen_db, int N_full,
    cudaStream_t stream,
    __half2* d_precomputed_query_rfft)  // Changed to use pre-allocated tensor
{
    // Allocate memory for precomputed query RFFT coefficients
    //__half2* d_precomputed_query_rfft;
    //int H = N_full >> 1;

    //printf("N_full: %d\n", N_full); crambin has N_full=64
    
    switch (N_full) {
        case 64: {
            int threads_per_block = 32; // 64 >> 1
            int query_shared_mem_size = 32 * sizeof(__half2); // H elements for query preprocessing
            //int main_shared_mem_size = (33 + 32 + 33 + 33) * sizeof(__half2); // (H+1) + H + (H+1) + (H+1)
            //int main_shared_mem_size = (33 + 32 + 33 + 33) * sizeof(__half2) + 3 * 64 * sizeof(float);
            int main_shared_mem_size = (33) * sizeof(__half2);
            
            // Precompute query RFFTs
            precompute_query_rfft_kernel<6><<<k, threads_per_block, query_shared_mem_size, stream>>>(
                d_precomputed_query_rfft,
                query_reversed_indices_ptr,
                weight_LUTS_A_ptr,
                k, vocab_size, query_len
            );
            
            // Main correlation kernel
            corr_fft_warp<6><<<B, threads_per_block, main_shared_mem_size, stream>>>(
                corr_ptr,
                d_precomputed_query_rfft,
                W_db_ptr,
                flat_db_ptr,
                offsets_db_ptr,
                lengths_db_ptr,
                k, vocab_size, maxlen_db
            );
            
            break;
        }
        case 128: {
            int threads_per_block = 64; // 128 >> 1
            int query_shared_mem_size = 64 * sizeof(__half2);
            // only use shared memory now to do the weird permutation thing. 
            int main_shared_mem_size = (65) * sizeof(__half2); 
            
            precompute_query_rfft_kernel<7><<<k, threads_per_block, query_shared_mem_size, stream>>>(
                d_precomputed_query_rfft,
                query_reversed_indices_ptr,
                weight_LUTS_A_ptr,
                k, vocab_size, query_len
            );

            threads_per_block = 32;//32;
            
            corr_fft_warp<7><<<B, threads_per_block, main_shared_mem_size, stream>>>(
                corr_ptr,
                d_precomputed_query_rfft,
                W_db_ptr,
                flat_db_ptr,
                offsets_db_ptr,
                lengths_db_ptr,
                k, vocab_size, maxlen_db
            );
            break;
        }
        case 256: {
            int threads_per_block = 128; // 256 >> 1
            int query_shared_mem_size = 128 * sizeof(__half2);
            //int main_shared_mem_size = (129 + 128 + 129 + 129) * sizeof(__half2);
            int main_shared_mem_size = (129) * sizeof(__half2);
            
            precompute_query_rfft_kernel<8><<<k, threads_per_block, query_shared_mem_size, stream>>>(
                d_precomputed_query_rfft,
                query_reversed_indices_ptr,
                weight_LUTS_A_ptr,
                k, vocab_size, query_len
            );
            threads_per_block = 32;
           
            corr_fft_warp<8><<<B, threads_per_block, main_shared_mem_size, stream>>>(
                corr_ptr,
                d_precomputed_query_rfft,
                W_db_ptr,
                flat_db_ptr,
                offsets_db_ptr,
                lengths_db_ptr,
                k, vocab_size, maxlen_db
            );
            break;
        }
        case 512: {
            int threads_per_block = 256; // 512 >> 1
            int query_shared_mem_size = 256 * sizeof(__half2);
            //int main_shared_mem_size = (257 + 256 + 257 + 257) * sizeof(__half2);
            int main_shared_mem_size = (257) * sizeof(__half2);
            
            precompute_query_rfft_kernel<9><<<k, threads_per_block, query_shared_mem_size, stream>>>(
                d_precomputed_query_rfft,
                query_reversed_indices_ptr,
                weight_LUTS_A_ptr,
                k, vocab_size, query_len
            );
            threads_per_block = 32;
          
            //corr_fft_warp<9><<<B, threads_per_block, 1, stream>>>( //@alex: 1->main_shared_mem_size appears to be bug
            corr_fft_warp<9><<<B, threads_per_block, main_shared_mem_size, stream>>>(
                corr_ptr,
                d_precomputed_query_rfft,
                W_db_ptr,
                flat_db_ptr,
                offsets_db_ptr,
                lengths_db_ptr,
                k, vocab_size, maxlen_db
            );
            break;
        }
        case 1024: {
            int threads_per_block = 512; // 1024 >> 1
            int query_shared_mem_size = 512 * sizeof(__half2);
            //int main_shared_mem_size = (513 + 512 + 513 + 513) * sizeof(__half2);
            int main_shared_mem_size = (513) * sizeof(__half2);
            
            precompute_query_rfft_kernel<10><<<k, threads_per_block, query_shared_mem_size, stream>>>(
                d_precomputed_query_rfft,
                query_reversed_indices_ptr,
                weight_LUTS_A_ptr,
                k, vocab_size, query_len
            );
            threads_per_block = 32;
            
            corr_fft_warp<10><<<B, threads_per_block, main_shared_mem_size, stream>>>(
                corr_ptr,
                d_precomputed_query_rfft,
                W_db_ptr,
                flat_db_ptr,
                offsets_db_ptr,
                lengths_db_ptr,
                k, vocab_size, maxlen_db
            );
            break;
        }
        case 2048: {
            int threads_per_block = 1024; // 1024 >> 1
            int query_shared_mem_size = 1024 * sizeof(__half2);
            //int main_shared_mem_size = (513 + 512 + 513 + 513) * sizeof(__half2);
            int main_shared_mem_size = (1025) * sizeof(__half2);

            //printf("2048\n");
            break;
            
            precompute_query_rfft_kernel<11><<<k, threads_per_block, query_shared_mem_size, stream>>>(
                d_precomputed_query_rfft,
                query_reversed_indices_ptr,
                weight_LUTS_A_ptr,
                k, vocab_size, query_len
            );
            printf("2048\n");
            threads_per_block = 32;
            
            corr_fft_warp<11><<<B, threads_per_block, main_shared_mem_size, stream>>>(
                corr_ptr,
                d_precomputed_query_rfft,
                W_db_ptr,
                flat_db_ptr,
                offsets_db_ptr,
                lengths_db_ptr,
                k, vocab_size, maxlen_db
            );
            printf("2048\n");
            break;
        }
        case 4096: {
            int threads_per_block = 2048; // 1024 >> 1
            int query_shared_mem_size = 2048 * sizeof(__half2);
            //int main_shared_mem_size = (513 + 512 + 513 + 513) * sizeof(__half2);
            int main_shared_mem_size = (2049) * sizeof(__half2);
            
            precompute_query_rfft_kernel<12><<<k, threads_per_block, query_shared_mem_size, stream>>>(
                d_precomputed_query_rfft,
                query_reversed_indices_ptr,
                weight_LUTS_A_ptr,
                k, vocab_size, query_len
            );
            threads_per_block = 32;
            
            corr_fft_warp<12><<<B, threads_per_block, main_shared_mem_size, stream>>>(
                corr_ptr,
                d_precomputed_query_rfft,
                W_db_ptr,
                flat_db_ptr,
                offsets_db_ptr,
                lengths_db_ptr,
                k, vocab_size, maxlen_db
            );
            break;
        }
        case 8192: {
            int threads_per_block = 4096; // 1024 >> 1
            int query_shared_mem_size = 4096 * sizeof(__half2);
            //int main_shared_mem_size = (513 + 512 + 513 + 513) * sizeof(__half2);
            int main_shared_mem_size = (4097) * sizeof(__half2);
            
            precompute_query_rfft_kernel<12><<<k, threads_per_block, query_shared_mem_size, stream>>>(
                d_precomputed_query_rfft,
                query_reversed_indices_ptr,
                weight_LUTS_A_ptr,
                k, vocab_size, query_len
            );
            threads_per_block = 32;
            
            corr_fft_warp<12><<<B, threads_per_block, main_shared_mem_size, stream>>>(
                corr_ptr,
                d_precomputed_query_rfft,
                W_db_ptr,
                flat_db_ptr,
                offsets_db_ptr,
                lengths_db_ptr,
                k, vocab_size, maxlen_db
            );
            break;
        }
    }
    
}