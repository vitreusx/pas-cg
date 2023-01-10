#pragma once
#include <cg/vect/vect.h>
#include <cuda.h>
#include <type_traits>
#include <utility>

namespace cg::cuda {
template <typename T, typename CPUAlloc = vect::allocator<T>,
          typename GPUAlloc = vect::cuda_allocator<T>>
class dev_var {
public:
  template <typename... Args>
  dev_var(Args &&...args) {
    *this = T(std::forward<Args>(args)...);
  }

  dev_var(dev_var const &other) {
    cpu_p = cpu_alloc.allocate(1);
    gpu_p = gpu_alloc.allocate(1);
    *this = static_cast<T>(other);
  }

  dev_var(dev_var &&other) {
    cpu_p = other.cpu_p;
    other.cpu_p = nullptr;
    gpu_p = other.gpu_p;
    other.gpu_p = nullptr;
  }

  ~dev_var() {
    destroy();
  }

  dev_var &operator=(dev_var const &other) {
    *this = static_cast<T>(other);
    return *this;
  }

  dev_var &operator=(dev_var &&other) {
    *this = static_cast<T>(other);
    other.destroy();
    return *this;
  }

  template <typename U>
  dev_var &operator=(U const &val) {
    set<U>(val);
    return *this;
  }

  template <typename U>
  void set(U const &val, cudaStream_t stream = 0) {
    if (!cpu_p)
      cpu_p = cpu_alloc.allocate(1);
    if (!gpu_p)
      gpu_p = gpu_alloc.allocate(1);

    *cpu_p = val;
    cudaMemcpyAsync(gpu_p, cpu_p, sizeof(T), cudaMemcpyHostToDevice, stream);
  }

  void fetch(cudaStream_t stream = 0) {
    cudaMemcpyAsync(cpu_p, gpu_p, sizeof(T), cudaMemcpyDeviceToHost, stream);
  }

  operator T() const {
    return *cpu_p;
  }

  T *get() {
    return gpu_p;
  }

  T const *get() const {
    return gpu_p;
  }

private:
  T *cpu_p = nullptr, *gpu_p = nullptr;
  CPUAlloc cpu_alloc;
  GPUAlloc gpu_alloc;

  void destroy() {
    if (cpu_p) {
      cpu_alloc.deallocate(cpu_p, 1);
      cpu_p = nullptr;
    }

    if (gpu_p) {
      gpu_alloc.deallocate(gpu_p, 1);
      gpu_p = nullptr;
    }
  }
};
} // namespace cg::cuda