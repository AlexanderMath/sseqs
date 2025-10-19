#include <torch/extension.h>
#include <ATen/cuda/CUDAContext.h> 
#include <cuda_runtime.h> 
#include <stdint.h> 

// Forward declarations
extern "C"
void launch_sw_cuda_blosum62(
        const uint8_t* query_seq_indices_ptr,
        int            query_length,
        int*           good_idx,
        const uint8_t* ascii,
        int*           lengths,
        int            num_sequences_in_db,
        int*           output_scores_ptr,
        int            gap_val,
        cudaStream_t   stream);
extern "C"
void launch_sw_cuda_affine(
        const uint8_t* query_seq_indices_ptr,
        int            query_length,
        int*           good_idx,
        const uint8_t* ascii,
        int*           lengths,
        int            num_sequences_in_db,
        int*           output_scores_ptr,
        int            gap_open,
        int            gap_extend,
        cudaStream_t   stream);


torch::Tensor sw_cuda_blosum62_pybind_wrapper(
        torch::Tensor query_indices_tensor,
        torch::Tensor good_idx,
        torch::Tensor ascii,
        torch::Tensor lengths,
        int           gap_penalty_param = -11) {
    // ---- Prepare Data for CUDA Kernel ----
    const int query_length = query_indices_tensor.size(0);
    const int num_db_sequences = good_idx.size(0); // offsets has nseq+1 elements

    // Allocate output tensor for scores
    auto output_scores_tensor = torch::empty({num_db_sequences},
                                             torch::TensorOptions()
                                                 .dtype(torch::kInt32)
                                                 .device(query_indices_tensor.device()));

    // Get raw data pointers
    const uint8_t* query_ptr   = query_indices_tensor.data_ptr<uint8_t>();
    int*           scores_ptr  = output_scores_tensor.data_ptr<int>();

    // Get current CUDA stream from PyTorch
     auto stream = at::cuda::getCurrentCUDAStream();
    //cudaStream_t stream = c10::cuda::getCurrentCUDAStream(query_indices_tensor.device().index());

    // ---- Launch CUDA Kernel ----
    launch_sw_cuda_blosum62(
        query_ptr,
        query_length,
        good_idx.data_ptr<int>(),
        ascii.data_ptr<uint8_t>(),
        lengths.data_ptr<int>(),
        num_db_sequences,
        scores_ptr,
        gap_penalty_param,
        stream
    );
    
    // cudaError_t err = cudaGetLastError(); // Optional: Check for errors after kernel launch
    // TORCH_CHECK(err == cudaSuccess, "CUDA error after sw_kernel launch: ", cudaGetErrorString(err));

    return output_scores_tensor;
}

// -----------------------------------------------------------------------------
//  PyTorch C++ wrapper function for Affine Smith-Waterman
// -----------------------------------------------------------------------------
torch::Tensor sw_cuda_affine_pybind_wrapper(
        torch::Tensor query_indices_tensor,
        torch::Tensor good_idx,
        torch::Tensor ascii,
        torch::Tensor lengths,
        int           gap_open = 11,
        int           gap_extend = 1) // Default gap penalties
{
    const int query_length = query_indices_tensor.size(0);
    const int num_db_sequences = good_idx.size(0);
    auto output_scores_tensor = torch::empty({num_db_sequences},
                                             torch::TensorOptions()
                                                 .dtype(torch::kInt32)
                                                 .device(query_indices_tensor.device()));
    const uint8_t* query_ptr   = query_indices_tensor.data_ptr<uint8_t>();
    int*           scores_ptr  = output_scores_tensor.data_ptr<int>();
    auto stream = at::cuda::getCurrentCUDAStream();
    launch_sw_cuda_affine(
        query_ptr,
        query_length,
        good_idx.data_ptr<int>(),
        ascii.data_ptr<uint8_t>(),
        lengths.data_ptr<int>(),
        num_db_sequences,
        scores_ptr,
        gap_open,
        gap_extend,
        stream
    );
    return output_scores_tensor;
}

// =============================== PYBIND11 MODULE DEFINITION ================================ //
PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) { 
    m.def("sw_cuda_blosum62", &sw_cuda_blosum62_pybind_wrapper,
            R"doc(Smith-Waterman (CUDA) with BLOSUM62 scoring.
            Returns SW score for each DB sequence as torch.Tensor (int32, 1D, CUDA))doc",
          py::arg("query_indices_tensor"),
          py::arg("good_idx"),
          py::arg("ascii_tensor"),
          py::arg("lengths_tensor"),
          py::arg("gap_penalty_param")      = -10
    );

    m.def("sw_cuda_affine", &sw_cuda_affine_pybind_wrapper,
          R"doc(Smith-Waterman (CUDA) with BLOSUM62 scoring and affine gap penalties.
          Returns affine SW score for each DB sequence as torch.Tensor (int32, 1D, CUDA))doc",
          py::arg("query_indices_tensor"),
          py::arg("good_idx"),
          py::arg("ascii_tensor"),
          py::arg("lengths_tensor"),
          py::arg("gap_open")      = 11,
          py::arg("gap_extend")    = 1
    );
}