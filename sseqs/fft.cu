#include <cuda_runtime.h>
#include "helper.cu"

// ── AA-to-index: 0-19 = "ARNDCQEGHILKMFPSTWYV", 20 = X/-/*, 21 = other ───────────
__device__ __constant__ unsigned char ASCII2POS[256] = {
    /* default 21 for all 256 slots, then overwrite the known codes */
    #define XX 21
    [0 ... 255] = XX,
    ['A'] = 0,  ['R'] = 1,  ['N'] = 2,  ['D'] = 3,
    ['C'] = 4,  ['Q'] = 5,  ['E'] = 6,  ['G'] = 7,
    ['H'] = 8,  ['I'] = 9,  ['L'] = 10, ['K'] = 11,
    ['M'] = 12, ['F'] = 13, ['P'] = 14, ['S'] = 15,
    ['T'] = 16, ['W'] = 17, ['Y'] = 18, ['V'] = 19,
    ['X'] = 20, ['-'] = 20, ['*'] = 20
};

template<int LOG_N>
__global__ void corr_fft_warp(
          __half*            __restrict__ corr,   // (B , N)
    const __half2*           __restrict__ precomputed_query_rfft, // (k, H+1) precomputed query RFFTs
    const float*           __restrict__ W_db,                 // (R_lut, vocab_size_padded) float values for database (_weight_LUTS_B)
    const uint8_t*         __restrict__ flat_db,              // DB sequence token indices
    const int32_t*         __restrict__ offsets_db,
    const int32_t*         __restrict__ lengths_db,
    int  k,            // rank (row in weight_LUTS_A_flat and W_db)
    int  vocab_size,   // vocab_size_padded (columns in LUTS_A and W_db)
    int  maxlen_db)    // clamp on lengths for DB sequences
{

    constexpr int  N      = 1 << LOG_N;       // time-domain length (power-of-2)
    constexpr int  H      = N >> 1;           // complex length after pairing (N/2)
    constexpr int  LOG_H  = LOG_N - 1;
    //constexpr float invN  = 1.0f / N;
    constexpr int  TPT    = H / 32;       // how many each lane owns
    const float M_PI_f = __uint_as_float(0x40490FDB);

    // precompute twiddles 
    __shared__ __half2 s_twiddle[2*H];          // size: 2*H
    __half2* s_Wstage = s_twiddle;                  // 0 … H-2
    __half2* s_WNk    = s_twiddle + (H - 1);        // H values
    const int lane = threadIdx.x;             // 0 … 31   (one warp / block)
    int off = 0;                              // running offset inside s_Wstage[]
    #pragma unroll
    for (int stage = 0; stage < LOG_H; ++stage) {
        const int half         = 1 << stage;              // 1,2,4,…
        const float k2pi_over_m = -2.0f * M_PI_f / (half << 1);

        // each lane writes j = lane, lane+32, lane+64, …  < half
        for (int j = lane; j < half; j += 32) {
            float θ = k2pi_over_m * j;
            s_Wstage[off + j] = h2(cosf(θ), sinf(θ));
        }
        off += half;                                      // advance to next block
    }
    for (int k_idx = lane; k_idx < H; k_idx += 32) {
        float θ = 2.0f * M_PI_f * k_idx / N;              // ± sign later via conj
        s_WNk[k_idx] = h2(cosf(θ), sinf(θ));
    }
    __syncwarp();      // twiddles precomputed. 


    /* -------------------------------------------------------------- */
    /* block/thread bookkeeping                                       */
    /* -------------------------------------------------------------- */
    int bid = blockIdx.x;        // sequence id (for DB)
    int tid = threadIdx.x;       // 0 … H-1  (we launch H threads)

    if (tid >= H) return;        // only H threads participate

    // Shared memory layout:
    __shared__ __half2 s_db_seq_packed_time_Gk_ifft_input[H+1];

    __half2 v[TPT];
    __half2 v_prod[TPT];  // 
    __half2 q[TPT];       
    __half2 accumulate[TPT] = {h2(0.f, 0.f)};
    __half2 qH;           
    uint32_t mask = __activemask();          // safer than a hard-coded constant
    constexpr unsigned ALL = 0xffffffffu;
    const int base = lane * TPT;          // first index this lane owns
    //float sum  = 0.f;
    const unsigned ACTIVE_MASK = __activemask();   
    const int lanes      = H / TPT;                
    const int lane_pairs = lanes - 1;              
    __half2 v0; 
    int seq_len_db = lengths_db[bid];
    int start_offset_db = offsets_db[bid];
    //int start_offset = start_offset_db;

    // could reduce bandwidth 2x with float->__half. 
    __shared__ float s_W_db[5*20]; // weirdly <5 breaks, confusing. 
    #pragma unroll
    for (int i = threadIdx.x; i < 5*20; i+=blockDim.x) s_W_db[i] = W_db[i]; 

    // load the sequence into shared memory. 
    __shared__ uint8_t s_flat_db[N];
    #pragma unroll
    for (int i = threadIdx.x; i < N; i+=blockDim.x) {
        if (i < seq_len_db) { 
            s_flat_db[i] = ASCII2POS[flat_db[start_offset_db+i]];
            //s_flat_db[i] = flat_db[start_offset_db+i];
        }
        else s_flat_db[i] = 0;
    }

    #pragma unroll
    for (int r = 0; r < k; r++) {
        // load query[r] to local variable q[].             0.22ms, loading q[j]=h(.0f,.0f) makes ~7.13->6.91ms 
        const __half2* query_rfft_row = precomputed_query_rfft + r * (H + 1);
        for (int j = 0; j < TPT; ++j) q[j] = query_rfft_row[tid*TPT + j];
        if (tid == 0)  qH = query_rfft_row[H];

        // load db_sequences[bid] to local variable v[]     (bid=blockIdx.x)
        // loading bit-reversed. 
        float w0, w1; 
        #pragma unroll
        for (int j = 0; j < TPT; ++j) {
            int g     = tid * TPT + j;                         // global index
            int brevd = __brev(g) >> (32 - LOG_H);             // bit-reversed
            int idx0  = 2 * brevd;                             // even token
           
            // adding 0 here makes it 8.83ms -> 7.7ms
            // after hcmadd, it goes 7.13ms -> 6.06ms    (-> 6.58 if we do w0=s_W_db[g] so ~0.5ms)
            if (idx0  < seq_len_db && idx0  < maxlen_db) w0 = s_W_db[r * vocab_size + s_flat_db[idx0]];
            else w0 = 0.f;
            if (idx0+1 < seq_len_db && idx0+1 < maxlen_db && idx0+1 < N) w1 = s_W_db[r * vocab_size + s_flat_db[idx0+1]];
            else w1 = 0.f;
            v[j] = h2(w0, w1);
        }

        // 2.6ms 
        //for (int j = 0; j < TPT; ++j) { accumulate[j] = __hadd2(accumulate[j], v[j]); }continue; 
 
        // code C&P from fft_inplace_warp see filter.cu 
        #pragma unroll
        for (int stage = 0; stage < LOG_H; ++stage) {
            const int m    = 1 << (stage + 1);   // 2,4,8,…
            const int half = m >> 1;             // 1,2,4,…
            int offset = (1<<stage)-1;

            if (half < TPT) {
                /* partner is another element **inside the same lane** */
                for (int i = 0; i < TPT; ++i) {
                    int g = base + i;
                    if ((g & half) == 0) {
                        int j = g & (half - 1);
                        __half2 w = s_Wstage[offset+j];
                        __half2 a = v[i];
                        __half2 b = v[i ^ half];
                        __half2 t = cmul_fast(b, w);
                        v[i]         = __hadd2(a, t);
                        v[i ^ half]  = __hsub2(a, t);
                    }
                }
            }
            else {
                /* partner lives in ANOTHER lane → one shuffle          */
                const int lane_delta = half / TPT;
                for (int i = 0; i < TPT; ++i) {
                    int g = base + i;
                    bool lower = (g & half) == 0;
                    __half2 a = v[i];
                    __half2 b = __shfl_xor_sync(ALL, a, lane_delta);
                    int   j = g & (half - 1);
                    __half2 w = s_Wstage[offset+j];
                    //if (lower) v[i] = __hadd2(a, cmul_fast(b, w));      
                    if (lower) v[i] = __hcmadd(b, w, a);
                    //else       v[i] = __hsub2(b, cmul_fast(a, w));      
                    else       v[i] = __hcmadd(__hneg2(a), w, b); // 8.10ms -> 8.05ms 
                }
            }
        }

        // accumulate in registers. 
        // 3.3ms -- only 0.7ms
        //for (int j = 0; j < TPT; ++j) { accumulate[j] = __hadd2(accumulate[j], v[j]); } continue; 
        
        #pragma unroll
        for (int j = 0; j < TPT; ++j) {
            int k_idx = tid * TPT + j;                 // 0 … H-1  (for live lanes)

            /* ---------- locate partner (H-k) inside the warp ----------------- */
            int mate_lane, mate_j;
            if (j == 0) {                              // special row
                mate_j    = 0;                         // slot-0 mirrors slot-0
                mate_lane = (tid == 0) ? 0            // k==0 lane mirrors itself
                                    : lanes - tid;  // all other lanes: H-k
            } else {                                   // j = 1 … TPT-1
                mate_j    = TPT - j;                   // 1↔TPT-1, 2↔TPT-2, …
                mate_lane = lane_pairs - tid;          // 0↔15, 1↔14, …
            }

            __half2 Gmate = __shfl_sync(ACTIVE_MASK, v[mate_j], mate_lane);
            __half2 GHmk  = conjh(Gmate);           

            if (k_idx == 0) {
                __half2 P0 = cmul_fast(q[0], h2(real(v[0]) + imag(v[0]), 0.f));
                __half2 PH = h2(real(cmul_fast(qH, h2(real(v[0]) - imag(v[0]), 0.f))), 0.f);
                v[0] = h2(0.5f * (real(P0) + real(PH)), 0.5f * (real(P0) - real(PH)));
                v0 = v[0];
            }
            else {
                __half2 Gk    = v[j];
                __half2 WNk   = conjh(s_WNk[k_idx]);
                __half2 term_sum  = __hadd2(Gk, GHmk);
                __half2 term_diff = __hsub2(Gk, GHmk);
                //__half2 jWNk_diff = cmul_fast(h2(0.f, -1.f), cmul_fast(WNk, term_diff));
                //v_prod[j] = cmul_fast(q[j], h2(0.5f * real(add(term_sum, jWNk_diff)), 0.5f * imag(add(term_sum, jWNk_diff))));

                // 7.39ms -> 7.28ms
                __half2 tmp = __hcmadd(h2(0.f, -1.f),           // -i
                       cmul_fast(WNk, term_diff), term_sum);               // a·b + c
                // 2)   v_prod[j] = q[j] · 0.5·tmp     (scale once, reuse both parts)
                v_prod[j] = cmul_fast(q[j], __hmul2(tmp, h2(0.5f, 0.5f)));
            }
        }

        // 3.79ms
        //for (int j = 0; j < TPT; ++j) { accumulate[j] = __hadd2(accumulate[j], v_prod[j]); } continue; 

        #pragma unroll
        for (int j = 0; j < TPT; ++j) {
            int k_idx = tid * TPT + j;               // logical k = 0 … H-1

            // partner location ---------------------------------------------------- 
            int mate_lane, mate_j;
            if (j == 0) {
                mate_j    = 0;
                mate_lane = (tid == 0) ? 0 : (lanes - tid);
            } else {
                mate_j    = TPT - j;
                mate_lane = lane_pairs - tid;
            }

            // every lane participates in the shuffle ----------------------------- 
            __half2 Pk    = v_prod[j];
            __half2 Pmate = __shfl_sync(mask, v_prod[mate_j], mate_lane);
            __half2 PHmk_c = conjh(Pmate);

            // skip the math/store for k = 0 (after the shuffle, so no deadlock) 
            if (k_idx == 0)       continue;

            // split-radix butterfly ---------------------------------------------- 
            __half2 WNnegk  = s_WNk[k_idx];

            __half2 term_sum  = __hadd2(Pk, PHmk_c);
            __half2 term_diff = __hsub2(Pk, PHmk_c);
            __half2 jWN_diff  = cmul_fast(h2(0.f, 1.f), cmul_fast(WNnegk, term_diff));
            __half2 tmp = h2(0.5f * real(add(term_sum, jWN_diff)), 0.5f * imag(add(term_sum, jWN_diff)));
            // wasn't faster with hcmadd
            //__half2 tmp  = __hcmadd(h2(0.f, 1.f),                           //  i
            //            cmul_fast(WNnegk, term_diff),           //  WNnegk·term_diff
            //            term_sum);                              //  + term_sum
            //tmp = __hmul2(tmp, h2(0.5f, 0.5f));                             // scale by ½

            s_db_seq_packed_time_Gk_ifft_input[k_idx] = tmp; // 0.5ms
            v[j] = tmp; 
        }

        // 4.73ms -- above is 1ms!   4.21ms if we remove the s_db_seq[k_idx]=tmp;
        //for (int j = 0; j < TPT; ++j) { accumulate[j] = __hadd2(accumulate[j], v[j]); } continue; 

        #pragma unroll
        for (int j = 0; j < TPT; ++j) {
            int rev         = __brev(tid*TPT+j) >> (32 - LOG_H);
            // swap each pair exactly once, just like the original code 
            if (tid*TPT+j < rev) {
                // using v[j] here instead of s_db_seq makes it 9.9ms -> 9ms. 
                // quite confused on interpreting the profile. 
                __half2 tmp = v[j];//s_db_seq_packed_time_Gk_ifft_input[tid*TPT+j];
                s_db_seq_packed_time_Gk_ifft_input[tid*TPT+j] = s_db_seq_packed_time_Gk_ifft_input[rev];
                s_db_seq_packed_time_Gk_ifft_input[rev] = tmp; // removing these do not speed up! 
            }
        }

        // 5.13ms -- above is 0.4ms! 
        //for (int j = 0; j < TPT; ++j) { accumulate[j] = __hadd2(accumulate[j], s_db_seq_packed_time_Gk_ifft_input[j]); } continue; 

        #pragma unroll
        for (int j = 0; j < TPT; ++j) {
            if (tid == 0 && j == 0) continue; 
            v[j] = s_db_seq_packed_time_Gk_ifft_input[tid*TPT + j];
        }

        // 5.22ms -- 0.1ms, fast read! 
        //for (int j = 0; j < TPT; ++j) { accumulate[j] = __hadd2(accumulate[j], v[j]); } continue; 

        
        // -------------- inv fft start -------------------
        #pragma unroll
        for (int stage = 0; stage < LOG_H; ++stage) {
            const int m    = 1 << (stage + 1);   // 2,4,8,…
            const int half = m >> 1;             // 1,2,4,…
            int offset = (1<<stage)-1;

            if (half < TPT) {
                /* partner is another element **inside the same lane** */
                for (int i = 0; i < TPT; ++i) {
                    int g = base + i;
                    if ((g & half) == 0) {
                        int j = g & (half - 1);
                        __half2 w = conjh(s_Wstage[offset+j]);

                        __half2 a = v[i];
                        __half2 b = v[i ^ half];
                        //__half2 t = cmul_fast(b, w);
                        //v[i]         = __hadd2(a, t);
                        //v[i ^ half]  = __hsub2(a, t); // 7.75ms -> 7.43ms -> 7.39ms by using one hcmadd instead of two!
                        __half2 vplus = __hcmadd(b, w, a);                 // a + b·w
                        v[i]        = vplus; 
                        v[i ^ half] = __hsub2(__hadd2(a, a), vplus); //__hcmadd(b, __hneg2(w), a);        // a – b·w   (because b·(–w) = –b·w)
                    }
                }
            }
            else {
                /* partner lives in ANOTHER lane → one shuffle          */
                const int lane_delta = half / TPT;
                for (int i = 0; i < TPT; ++i)
                {
                    int g = base + i;
                    bool lower = (g & half) == 0;

                    __half2 a = v[i];
                    __half2 b = __shfl_xor_sync(ALL, a, lane_delta);

                    int   j = g & (half - 1);
                    __half2 w = conjh(s_Wstage[offset+j]);

                    //if (lower) v[i] = __hadd2(a, cmul_fast(b, w));      
                    if (lower) v[i] = __hcmadd(b, w, a); // 8.23ms -> 8.10ms
                    //else       v[i] = __hsub2(b, cmul_fast(a, w));      
                    else       v[i] = __hcmadd(__hneg2(a), w, b); // 8.05ms -> 7.75ms
                }
            }
        }
        
        // -------------- inv fft end -------------------
        // accumulate in registers. 
        #pragma unroll 
        for (int j = 0; j < TPT; ++j) { accumulate[j] = __hadd2(accumulate[j], v[j]); }

        // 8.8ms -- 3.6ms for ifft. can read and do fwd fft in 1.3ms ?!
        // 8.35ms if we do #pragma unroll. 
    }


    //__half* row = corr + bid * N;
    __half* row = corr+bid; 
    /*{
        const float current_scale = invN;
        __half* row = corr + bid * N;                // pointer to this sequence's row

        #pragma unroll
        for (int j = 0; j < TPT; ++j) {
            const int logical_tid = tid * TPT + j;          // 0 … H-1
            __half2 vj = accumulate[j];//v[j];//s_db_seq_packed_time_Gk_ifft_input[logical_tid];

            // scale by 1/N -------------------------------------------------- 
            vj = h2(real(vj) * current_scale, imag(vj) * current_scale);

            // write two real samples --------------------------------------- 
            const float xr  = real(vj);
            const float xi  = imag(vj);

            const int out0 = 2 * logical_tid;    // 0,2,4,…,510
            const int out1 = out0 + 1;           // 1,3,5,…,511

            if (out0 < N) row[out0] = __float2half(__half2float(row[out0]) + xr);
            if (out1 < N) row[out1] = __float2half(__half2float(row[out1]) + xi);
        }
    }

    __syncthreads();*/

    /* pointer to this sequence's output row ----------------------------- */
    /*******************************************************************
    * 1.  Per-thread reduction over its own TPT samples  (real part)
    *******************************************************************/
    float local_max = -FLT_MAX;
    const float scale = 1.0f / static_cast<float>(N);

    #pragma unroll
    for (int j = 0; j < TPT; ++j)
    {
        __half2 v = accumulate[j];          // ← correct index
        v = __halves2half2(__hmul(__low2half(v),  __float2half(scale)),
                        __hmul(__high2half(v), __float2half(scale)));

        float xr = __half2float(__low2half(v));   // real part
        local_max = fmaxf(local_max, xr);
    }

    /*******************************************************************
    * 2.  Cross-thread reduction in shared memory
    *******************************************************************/
    extern __shared__ float sdata[];     // one float per thread (32)
    sdata[tid] = local_max;
    __syncthreads();

    /* binary-tree reduction: 32 → 16 → 8 → … → 1  (shared mem only)   */
    for (int stride = blockDim.x >> 1; stride > 0; stride >>= 1)
    {
        if (tid < stride)
            sdata[tid] = fmaxf(sdata[tid], sdata[tid + stride]);
        __syncthreads();
    }

    /*******************************************************************
    * 3.  Thread 0 writes the final result once to corr[bid, 0]
    *******************************************************************/
    if (tid == 0)
        row[0] = __float2half(sdata[0]);

   
}