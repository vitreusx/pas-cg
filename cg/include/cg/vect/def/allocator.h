#pragma once
#include "../utils/cuda.h"
#include <cstdint>
#include <cstdlib>
#include <cuda_runtime_api.h>
#include <memory>

namespace nitro::def {
template <typename T>
struct allocator {
  allocator() {
    if (cudaGetDeviceCount(&num_devices) != cudaSuccess)
      num_devices = 0;
  }

  T *allocate(std::size_t n, const void * = nullptr) {
    //    std::size_t alignment = 512 >> 3, nbytes = n * sizeof(T);
    //    nbytes = alignment * ((nbytes + alignment - 1) / alignment);
    //    return (T *)std::aligned_alloc(alignment, nbytes);

    std::size_t alignment = 512 >> 3;
    if (num_devices > 0) {
      std::size_t space = n * sizeof(T) + alignment;
      T *buf;
      cudaMallocHost((void **)&buf, space);
      return (T *)std::align(alignment, n * sizeof(T), (void *&)buf, space);
    } else {
      std::size_t space = n * sizeof(T);
      space = alignment * ((space + alignment - 1) / alignment);
      return (T *)std::aligned_alloc(alignment, space);
    }
  }

  void deallocate(T *p, std::size_t) {
    if (num_devices > 0) {
      cudaFreeHost(p);
    } else {
      free((void *)p);
    }
  }

private:
  int num_devices = 0;
};
} // namespace nitro::def