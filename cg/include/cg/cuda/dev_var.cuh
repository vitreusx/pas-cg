#pragma once
#include <cuda.h>
#include <type_traits>
#include <utility>

namespace cg::cuda {
template <typename T>
class dev_var {
public:
  template <typename... Args>
  dev_var(Args &&...args) {
    cudaMalloc(&ptr, sizeof(T));
    *this = T(std::forward<Args>(args)...);
  }

  template <typename U>
  dev_var &operator=(U const &val) {
    T value = static_cast<T>(val);
    cudaMemcpy(ptr, &value, sizeof(T), cudaMemcpyHostToDevice);
    return *this;
  }

  operator T() const {
    alignas(T) char value_buf[sizeof(T)];
    cudaMemcpy(value_buf, &ptr, sizeof(T), cudaMemcpyDeviceToHost);
    return *reinterpret_cast<T const *>(value_buf);
  }

private:
  T *ptr;
};
} // namespace cg::cuda