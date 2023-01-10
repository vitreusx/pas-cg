#pragma once
#include <cuda_runtime_api.h>
#include <optional>

namespace cg::cuda {
class stream {
public:
  stream() {
    cudaStreamCreate(&handle);
  }

  stream(stream const &) = delete;

  stream(stream &&other) {
    handle = other.handle;
    other.handle = 0;
  }

  ~stream() {
    destroy();
  }

  stream &operator=(stream const &) = delete;

  stream &operator=(stream &&other) {
    destroy();
    handle = other.handle;
    other.handle = 0;
    return *this;
  }

  operator cudaStream_t const &() const {
    return handle;
  }

  operator cudaStream_t &() {
    return handle;
  }

private:
  cudaStream_t handle;

  void destroy() {
    if (handle)
      cudaStreamDestroy(handle);
  }
};
} // namespace cg::cuda