#pragma once
#include <cuda_runtime_api.h>
#include <stdexcept>

#define cudaCheck(expr)                                                        \
  {                                                                            \
    cudaError_t ret = (expr);                                                  \
    if (ret != cudaSuccess)                                                    \
      throw std::runtime_error(cudaGetErrorString(ret));                       \
  }
