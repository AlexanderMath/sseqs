#include <torch/extension.h>
#include <ATen/cuda/CUDAContext.h> 
#include <cuda_runtime.h> 
#include <stdint.h> 

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
        int            db_seq_length,
        cudaStream_t   stream);

extern "C"
void launch_sw_cuda_affine_uint8(
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
        cudaStream_t   stream);
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
        cudaStream_t   stream);

extern "C"
void launch_sw_cuda_linear(
        const uint8_t* query_seq_indices_ptr,
        int            query_length,
        int*           good_idx,
        const uint8_t* ascii,
        int*           lengths,
        int            num_sequences_in_db,
        int*           output_scores_ptr,
        int            gap_penalty,
        int            db_seq_length,
        cudaStream_t   stream);

extern "C"
void launch_sw_cuda_profile(
        const int8_t*  pssm_ptr,
        int            query_length,
        int*           good_idx,
        const uint8_t* ascii,
        int*           lengths,
        int            num_sequences_in_db,
        int*           output_scores_ptr,
        int            gap_open,
        int            gap_extend,
        int            db_seq_length,
        cudaStream_t   stream);



// -----------------------------------------------------------------------------
//  PyTorch C++ wrapper function for Affine Smith-Waterman
// -----------------------------------------------------------------------------
torch::Tensor sw_cuda_affine_pybind_wrapper(
        torch::Tensor query_indices_tensor,
        torch::Tensor good_idx,
        torch::Tensor ascii,
        torch::Tensor lengths,
        int           db_seq_length,  // NEW: sequence length from Python
        int           gap_open = 11,
        int           gap_extend = 1) // Default gap penalties
{
    const int query_length = query_indices_tensor.size(0);
    const int num_db_sequences = good_idx.size(0);
    auto output_scores_tensor = torch::zeros({num_db_sequences},
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
        db_seq_length,  // Pass it through
        stream
    );
    
    return output_scores_tensor;
}

// -----------------------------------------------------------------------------
//  PyTorch C++ wrapper function for Affine Smith-Waterman (uint8 version)
// -----------------------------------------------------------------------------
torch::Tensor sw_cuda_affine_uint8_pybind_wrapper(
        torch::Tensor query_indices_tensor,
        torch::Tensor good_idx,
        torch::Tensor ascii,
        torch::Tensor lengths,
        int           gap_open = 11,
        int           gap_extend = 1)
{
    const int query_length = query_indices_tensor.size(0);
    const int num_db_sequences = good_idx.size(0);
    auto output_scores_tensor = torch::zeros({num_db_sequences},
                                             torch::TensorOptions()
                                                 .dtype(torch::kInt32)
                                                 .device(query_indices_tensor.device()));
    const uint8_t* query_ptr   = query_indices_tensor.data_ptr<uint8_t>();
    int*           scores_ptr  = output_scores_tensor.data_ptr<int>();
    auto stream = at::cuda::getCurrentCUDAStream();
    launch_sw_cuda_affine_uint8(
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

// -----------------------------------------------------------------------------
//  PyTorch C++ wrapper function for Affine Smith-Waterman WITH BACKTRACKING
// -----------------------------------------------------------------------------
std::tuple<torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor> sw_cuda_affine_backtrack_pybind_wrapper(
        torch::Tensor query_indices_tensor,
        torch::Tensor good_idx,
        torch::Tensor ascii,
        torch::Tensor lengths,
        int           max_align_len = 2048,
        int           gap_open = 11,
        int           gap_extend = 1)
{
    const int query_length = query_indices_tensor.size(0);
    const int num_db_sequences = good_idx.size(0);
    
    auto output_scores_tensor = torch::zeros({num_db_sequences},
                                             torch::TensorOptions()
                                                 .dtype(torch::kInt32)
                                                 .device(query_indices_tensor.device()));
    
    auto output_align_lens_tensor = torch::zeros({num_db_sequences},
                                                   torch::TensorOptions()
                                                       .dtype(torch::kInt32)
                                                       .device(query_indices_tensor.device()));
    
    auto output_end_i_tensor = torch::zeros({num_db_sequences},
                                             torch::TensorOptions()
                                                 .dtype(torch::kInt32)
                                                 .device(query_indices_tensor.device()));
    
    auto output_end_j_tensor = torch::zeros({num_db_sequences},
                                             torch::TensorOptions()
                                                 .dtype(torch::kInt32)
                                                 .device(query_indices_tensor.device()));
    
    auto output_align_ops_tensor = torch::zeros({num_db_sequences, max_align_len},
                                                 torch::TensorOptions()
                                                     .dtype(torch::kInt8)
                                                     .device(query_indices_tensor.device()));
    
    const uint8_t* query_ptr   = query_indices_tensor.data_ptr<uint8_t>();
    int*           scores_ptr  = output_scores_tensor.data_ptr<int>();
    int*           align_lens_ptr = output_align_lens_tensor.data_ptr<int>();
    int*           end_i_ptr = output_end_i_tensor.data_ptr<int>();
    int*           end_j_ptr = output_end_j_tensor.data_ptr<int>();
    char*          align_ops_ptr = (char*)output_align_ops_tensor.data_ptr<int8_t>();
    auto stream = at::cuda::getCurrentCUDAStream();
    
    launch_sw_cuda_affine_backtrack(
        query_ptr,
        query_length,
        good_idx.data_ptr<int>(),
        ascii.data_ptr<uint8_t>(),
        lengths.data_ptr<int>(),
        num_db_sequences,
        scores_ptr,
        align_lens_ptr,
        end_i_ptr,
        end_j_ptr,
        align_ops_ptr,
        max_align_len,
        gap_open,
        gap_extend,
        stream
    );
    
    return std::make_tuple(output_scores_tensor, output_align_lens_tensor, output_end_i_tensor, output_end_j_tensor, output_align_ops_tensor);
}

// -----------------------------------------------------------------------------
//  PyTorch C++ wrapper function for Anti-Diagonal Affine Smith-Waterman
// -----------------------------------------------------------------------------
torch::Tensor sw_cuda_affine_antidiag_pybind_wrapper(
        torch::Tensor query_indices_tensor,
        torch::Tensor good_idx,
        torch::Tensor ascii,
        torch::Tensor lengths,
        int           gap_open = 11,
        int           gap_extend = 1)
{
    const int query_length = query_indices_tensor.size(0);
    const int num_db_sequences = good_idx.size(0);
    auto output_scores_tensor = torch::zeros({num_db_sequences},
                                             torch::TensorOptions()
                                                 .dtype(torch::kInt32)
                                                 .device(query_indices_tensor.device()));
    const uint8_t* query_ptr   = query_indices_tensor.data_ptr<uint8_t>();
    int*           scores_ptr  = output_scores_tensor.data_ptr<int>();
    auto stream = at::cuda::getCurrentCUDAStream();
    launch_sw_cuda_affine_antidiag(
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

// -----------------------------------------------------------------------------
//  PyTorch C++ wrapper function for Linear Gap Smith-Waterman (screening)
// -----------------------------------------------------------------------------
torch::Tensor sw_cuda_linear_pybind_wrapper(
        torch::Tensor query_indices_tensor,
        torch::Tensor good_idx,
        torch::Tensor ascii,
        torch::Tensor lengths,
        int           db_seq_length,
        int           gap_penalty = 1)  // Default linear gap penalty
{
    const int query_length = query_indices_tensor.size(0);
    const int num_db_sequences = good_idx.size(0);
    auto output_scores_tensor = torch::zeros({num_db_sequences},
                                             torch::TensorOptions()
                                                 .dtype(torch::kInt32)
                                                 .device(query_indices_tensor.device()));
    const uint8_t* query_ptr   = query_indices_tensor.data_ptr<uint8_t>();
    int*           scores_ptr  = output_scores_tensor.data_ptr<int>();
    auto stream = at::cuda::getCurrentCUDAStream();
    launch_sw_cuda_linear(
        query_ptr,
        query_length,
        good_idx.data_ptr<int>(),
        ascii.data_ptr<uint8_t>(),
        lengths.data_ptr<int>(),
        num_db_sequences,
        scores_ptr,
        gap_penalty,
        db_seq_length,
        stream
    );
    return output_scores_tensor;
}

// -----------------------------------------------------------------------------
//  PyTorch C++ wrapper function for Profile-based Smith-Waterman
// -----------------------------------------------------------------------------
torch::Tensor sw_cuda_profile_pybind_wrapper(
        torch::Tensor pssm_tensor,      // PSSM: (query_len, 20) int8
        torch::Tensor good_idx,
        torch::Tensor ascii,
        torch::Tensor lengths,
        int           db_seq_length,
        int           gap_open = 11,
        int           gap_extend = 1)
{
    const int query_length = pssm_tensor.size(0);
    const int num_db_sequences = good_idx.size(0);
    auto output_scores_tensor = torch::zeros({num_db_sequences},
                                             torch::TensorOptions()
                                                 .dtype(torch::kInt32)
                                                 .device(pssm_tensor.device()));
    const int8_t* pssm_ptr    = pssm_tensor.data_ptr<int8_t>();
    int*          scores_ptr  = output_scores_tensor.data_ptr<int>();
    auto stream = at::cuda::getCurrentCUDAStream();
    launch_sw_cuda_profile(
        pssm_ptr,
        query_length,
        good_idx.data_ptr<int>(),
        ascii.data_ptr<uint8_t>(),
        lengths.data_ptr<int>(),
        num_db_sequences,
        scores_ptr,
        gap_open,
        gap_extend,
        db_seq_length,
        stream
    );
    return output_scores_tensor;
}

// =============================== PYBIND11 MODULE DEFINITION ================================ //
PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) { 

    m.def("sw_cuda_affine", &sw_cuda_affine_pybind_wrapper,
          R"doc(Smith-Waterman (CUDA) with BLOSUM62 scoring and affine gap penalties.
          Returns affine SW score for each DB sequence as torch.Tensor (int32, 1D, CUDA))doc",
          py::arg("query_indices_tensor"),
          py::arg("good_idx"),
          py::arg("ascii_tensor"),
          py::arg("lengths_tensor"),
          py::arg("db_seq_len"),
          py::arg("gap_open")      = 11,
          py::arg("gap_extend")    = 1
    );

    m.def("sw_cuda_affine_uint8", &sw_cuda_affine_uint8_pybind_wrapper,
          R"doc(Smith-Waterman (CUDA) with BLOSUM62 scoring and affine gap penalties (uint8 version for screening).
          Returns affine SW score for each DB sequence as torch.Tensor (int32, 1D, CUDA))doc",
          py::arg("query_indices_tensor"),
          py::arg("good_idx"),
          py::arg("ascii_tensor"),
          py::arg("lengths_tensor"),
          py::arg("gap_open")      = 11,
          py::arg("gap_extend")    = 1
    );

    m.def("sw_cuda_affine_backtrack", &sw_cuda_affine_backtrack_pybind_wrapper,
          R"doc(Smith-Waterman (CUDA) with BLOSUM62 scoring, affine gap penalties, and backtracking with checkpointing.)doc",
          py::arg("query_indices_tensor"),
          py::arg("good_idx"),
          py::arg("ascii_tensor"),
          py::arg("lengths_tensor"),
          py::arg("max_align_len") = 2048,
          py::arg("gap_open")      = 11,
          py::arg("gap_extend")    = 1
    );

    m.def("sw_cuda_affine_antidiag", &sw_cuda_affine_antidiag_pybind_wrapper,
          R"doc(Smith-Waterman (CUDA) with BLOSUM62 scoring and affine gap penalties, using anti-diagonal parallelization.
          This version processes anti-diagonals instead of rows, enabling parallel computation across each anti-diagonal.)doc",
          py::arg("query_indices_tensor"),
          py::arg("good_idx"),
          py::arg("ascii_tensor"),
          py::arg("lengths_tensor"),
          py::arg("gap_open")      = 11,
          py::arg("gap_extend")    = 1
    );

    m.def("sw_cuda_linear", &sw_cuda_linear_pybind_wrapper,
          R"doc(Smith-Waterman (CUDA) with BLOSUM62 scoring and linear gap penalties (for screening/prefilter).
          Returns linear SW score for each DB sequence as torch.Tensor (int32, 1D, CUDA))doc",
          py::arg("query_indices_tensor"),
          py::arg("good_idx"),
          py::arg("ascii_tensor"),
          py::arg("lengths_tensor"),
          py::arg("db_seq_len"),
          py::arg("gap_penalty")   = 1
    );

    m.def("sw_cuda_profile", &sw_cuda_profile_pybind_wrapper,
          R"doc(Smith-Waterman (CUDA) with PSSM (profile) scoring and affine gap penalties.
          Uses position-specific scoring matrix instead of BLOSUM62.
          PSSM shape: (query_len, 20) where each row contains scores for 20 amino acids.
          Returns profile SW score for each DB sequence as torch.Tensor (int32, 1D, CUDA))doc",
          py::arg("pssm_tensor"),
          py::arg("good_idx"),
          py::arg("ascii_tensor"),
          py::arg("lengths_tensor"),
          py::arg("db_seq_len"),
          py::arg("gap_open")      = 11,
          py::arg("gap_extend")    = 1
    );
}