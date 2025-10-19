#include <cuda.h>
#include <cuda_runtime.h>
#include <stdint.h>
#include <cstdio>    

// @alex: may be faster moving into fft kernel, would skip a VRAM<->registers copy.
// nvprof says it gets 300GB/s vs 24GB/s RAM->VRAM so 0.16s vs 3.6s negligable. 
extern "C"
__global__ void decompress_kernel(const uint8_t* __restrict__ in,
                                   uint8_t*       __restrict__ out,
                                   int            n_syms,
                                   int            lo)        // multiple of 8
{
    constexpr int OFFS = 64;                         // add back after unpack

    int t = blockIdx.x * blockDim.x + threadIdx.x;   // thread index
    int idx8 = t * 8;                                // first symbol this thread writes
    if (idx8 >= n_syms) return;

    /* load 5 packed bytes (40 bits) ---------------------------- */
    const uint8_t* p = in + t * 5;
    uint64_t window = uint64_t(p[0]) << 32 |
                      uint64_t(p[1]) << 24 |
                      uint64_t(p[2]) << 16 |
                      uint64_t(p[3]) <<  8 |
                      uint64_t(p[4]);

    /* unpack 8 symbols ---------------------------------------- */
#pragma unroll
    for (int k = 0; k < 8; ++k) {
        uint8_t sym = (window >> (35 - 5*k)) & 0x1F;   // 5 bits
        out[idx8 + k] = sym + OFFS;                    // restore ASCII
    }
}

extern "C"
__global__ void compress_kernel(const uint8_t* __restrict__ in,
                                 uint8_t*       __restrict__ out,
                                 int            n_syms)        // multiple of 8
{
    constexpr int OFFS = 64;      // subtract before pack

    int t = blockIdx.x * blockDim.x + threadIdx.x;
    int idx8 = t * 8;
    if (idx8 >= n_syms) return;

    /* gather 8 symbols and build 40-bit window ---------------- */
    uint64_t window = 0;
#pragma unroll
    for (int k = 0; k < 8; ++k) {
        uint8_t sym = in[idx8 + k] - OFFS;     // 0…31
        window |= uint64_t(sym & 0x1F) << (35 - 5*k);
    }

    /* write 5 packed bytes ----------------------------------- */
    uint8_t* p = out + t * 5;
    p[0] = uint8_t(window >> 32);
    p[1] = uint8_t(window >> 24);
    p[2] = uint8_t(window >> 16);
    p[3] = uint8_t(window >>  8);
    p[4] = uint8_t(window);
}

//  C-linkage launchers  
extern "C"
void decompress_launcher(const uint8_t* d_in, uint8_t* d_out, int n_syms, int lo) {
    const int BLOCK = 256;
    dim3 grid((n_syms/8 + BLOCK - 1) / BLOCK);
    d_in = d_in + lo; // d_in[lo:]
    decompress_kernel<<<grid, BLOCK, 0>>>(d_in, d_out, n_syms, lo);
}

extern "C"
void compress_launcher(const uint8_t* d_in, uint8_t* d_out, int n_syms) {
    const int BLOCK = 256;
    dim3 grid((n_syms/8 + BLOCK - 1) / BLOCK);
    compress_kernel<<<grid, BLOCK, 0>>>(d_in, d_out, n_syms);
}