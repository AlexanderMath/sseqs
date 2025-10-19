#include <torch/extension.h>
#include <stdint.h>

// launchers implemented in xbit.cu
extern "C" void compress_launcher  (const uint8_t*, uint8_t*, int);
extern "C" void decompress_launcher(const uint8_t*, uint8_t*, int, int);

// -----------------------------------------------------------------------------
// compression wrapper - hardcoded to 5 bits
// -----------------------------------------------------------------------------
torch::Tensor compress_cuda(torch::Tensor input) {
    input = input.contiguous().to(torch::kCUDA, torch::kUInt8);
    const int n = input.size(0);
    const int out_bytes = (n * 5 + 7) / 8;
    auto output = torch::zeros({out_bytes}, torch::dtype(torch::kUInt8).device(input.device()));
    compress_launcher(input.data_ptr<uint8_t>(), output.data_ptr<uint8_t>(), n);
    return output;
}

// -----------------------------------------------------------------------------
// decompression wrapper - hardcoded to 5 bits and offset 64
// -----------------------------------------------------------------------------
torch::Tensor decompress_cuda(torch::Tensor input, int n, int lo) {
    input = input.to(torch::kCUDA, torch::kUInt8);
    auto output = torch::empty({n}, torch::dtype(torch::kUInt8).device(input.device()));
    decompress_launcher(input.data_ptr<uint8_t>(), output.data_ptr<uint8_t>(), n, lo);
    return output;
}

// -----------------------------------------------------------------------------
// pybind11 module
// -----------------------------------------------------------------------------
PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {
    m.def("compress",   &compress_cuda,   "xbit compress (CUDA)");
    m.def("decompress", &decompress_cuda, "xbit decompress (CUDA)");
}