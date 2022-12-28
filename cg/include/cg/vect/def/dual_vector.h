#pragma once
#include "allocator.h"
#include "cuda_allocator.h"
#include "vector.h"

namespace nitro::def {
template <typename T, typename Alloc = allocator<T>,
          typename CudaAlloc = cuda_allocator<T>>
class dual_vector : public vector<T, Alloc> {
public:
  dual_vector() = default;

  explicit dual_vector(Alloc alloc, CudaAlloc cuda_alloc = CudaAlloc())
      : vector<T, Alloc>(std::move(alloc)), cuda_alloc(std::move(cuda_alloc)) {}

  explicit dual_vector(int n, T const &init = T(), Alloc alloc = Alloc(),
                       CudaAlloc cuda_alloc = CudaAlloc())
      : vector<T, Alloc>(n, init, std::move(alloc)),
        cuda_alloc(std::move(cuda_alloc)) {}

  void push_to_gpu() {
#ifdef __CUDA__
    if (cuda_capacity < this->_size) {
      if (cuda_buf)
        cuda_alloc.deallocate(cuda_buf, cuda_capacity);
      cuda_buf = cuda_alloc.allocate(this->_size);
      cuda_capacity = this->_size;
    }
    cudaMemcpy(cuda_buf, this->data, this->_size * sizeof(T),
               cudaMemcpyHostToDevice);
#endif
  }

  void pull_from_gpu() {
#ifdef __CUDA__
    cudaMemcpy(this->data, cuda_buf, this->_size * sizeof(T),
               cudaMemcpyDeviceToHost);
#endif
  }

  view<T> gpu_view() {
    return view<T>(cuda_buf, this->_size);
  }

  const_view<T> gpu_view() const {
    return const_view<T>(cuda_buf, this->_size);
  }

protected:
  void destroy() override {
#ifdef __CUDA__
    if (cuda_buf) {
      cuda_alloc.deallocate(cuda_buf, cuda_capacity);
      cuda_buf = nullptr;
      cuda_capacity = 0;
    }
#endif
    vector<T, Alloc>::destroy();
  }

  T *cuda_buf = nullptr;
  int cuda_capacity = 0;
  CudaAlloc cuda_alloc;
};
} // namespace nitro::def