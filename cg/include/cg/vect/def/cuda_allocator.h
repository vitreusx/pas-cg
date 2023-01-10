#pragma once
#include <cstdint>
#include <cuda_runtime_api.h>

namespace nitro::def {
template <typename T>
struct cuda_allocator {
  __host__ T *allocate(std::size_t n, const void * = nullptr) {
    T *buf;
    cudaMalloc((void **)&buf, n * sizeof(T));
    return buf;
  }

  __host__ void deallocate(T *p, std::size_t) {
    cudaFree(p);
  }
};
} // namespace nitro::def