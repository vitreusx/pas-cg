#pragma once
#include "../config.h"
#include <cstdint>

namespace nitro::def {
template <typename T>
struct cuda_allocator {
  __host__ T *allocate(std::size_t n, const void * = nullptr) {
#ifdef __CUDA__
    T *buf;
    cudaMallocHost(&buf, n * sizeof(T));
    return buf;
#endif
  }

  __host__ void deallocate(T *p, std::size_t) {
#ifdef __CUDA__
    cudaFreeHost(p);
#endif
  }
};
} // namespace nitro::def