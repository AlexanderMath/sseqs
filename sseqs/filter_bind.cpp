#include <torch/extension.h>
#include <ATen/cuda/CUDAContext.h>
#include <cuda_fp16.h> 

extern "C"
void launch_warp_corr(
    __half* corr_ptr,
    const uint8_t* query_reversed_indices_ptr,
    const float*   weight_LUTS_A_ptr,
    const float*   W_ptr,
    const uint8_t* flat_ptr,
    const int32_t* offsets_ptr,
    const int32_t* lengths_ptr,
    float*         out_stats_ptr,
    int B, int vocab_size, int k, int query_len, 
    int maxlen, int N_full,
    cudaStream_t stream,
    __half2* query_rfft_tensor_ptr);

torch::Tensor warp_corr(
                     torch::Tensor corr,    // (B,N) float32 or float16
                     torch::Tensor query_reversed_indices,
                     torch::Tensor weight_LUTS_A,
                     torch::Tensor W,       // (R,vocab) float32
                     torch::Tensor flat,    // uint8
                     torch::Tensor offsets, // int32 (B+1)
                     torch::Tensor lengths, // int32 (B)
                     int64_t        k,
                     int64_t        query_length,
                     int64_t        maxlen,
                     int64_t        n,
                     torch::Tensor query_rfft_tensor) {
    // Convert to half if it's float32
    torch::Tensor corr_half = corr;
    if (corr.scalar_type() == torch::kFloat32) corr_half = corr.to(torch::kFloat16);
    int N_full = n; 

    int B          = static_cast<int>(lengths.size(0));
    int vocab_size = static_cast<int>(W.size(1)); // vocab_size_padded
    int q_len      = static_cast<int>(query_reversed_indices.size(0));

    auto opts = torch::TensorOptions().dtype(torch::kFloat32).device(corr.device());
    auto out  = torch::empty({B, 2}, opts); 

    auto stream = at::cuda::getCurrentCUDAStream();

    launch_warp_corr(
        reinterpret_cast<__half*>(corr_half.data_ptr<at::Half>()),
        query_reversed_indices.data_ptr<uint8_t>(),
        weight_LUTS_A.data_ptr<float>(),
        W.data_ptr<float>(),
        flat.data_ptr<uint8_t>(),
        offsets.data_ptr<int32_t>(),
        lengths.data_ptr<int32_t>(),
        out.data_ptr<float>(), 
        B, vocab_size, static_cast<int>(k),
        q_len, // query_len for the kernel
        static_cast<int>(maxlen),
        N_full,
        stream,
        reinterpret_cast<__half2*>(query_rfft_tensor.data_ptr<at::Half>()));

    return out;
}

/* pybind */
PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {
    m.def("warp_corr", &warp_corr,
          "corr RFFT, queries FFT internal using a single warp per sequence. ",
          pybind11::arg("corr"),
          pybind11::arg("query_reversed_indices"),
          pybind11::arg("weight_LUTS_A"),
          pybind11::arg("W"),
          pybind11::arg("flat"),
          pybind11::arg("offsets"), 
          pybind11::arg("lengths"),
          pybind11::arg("k"), 
          pybind11::arg("query_length"), 
          pybind11::arg("maxlen"),
          pybind11::arg("n"),
          pybind11::arg("query_rfft_tensor"));
}
