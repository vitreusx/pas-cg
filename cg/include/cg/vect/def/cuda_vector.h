#pragma once
#include "cuda_allocator.h"
#include "vector.h"

namespace nitro::def {
template <typename T, typename Alloc = cuda_allocator<T>>
class cuda_vector {
public:
  cuda_vector() : data_{nullptr}, _size{0}, _capacity{0}, alloc{Alloc()} {}

  explicit cuda_vector(Alloc alloc)
      : data_{nullptr}, _size{0}, _capacity{0}, alloc{std::move(alloc)} {}

  template <typename CPUAlloc = allocator<T>>
  void pull_from(vector<T, CPUAlloc> const &v, cudaStream_t stream = 0) {
    reserve(v.size());
    cudaMemcpyAsync(data_, v.data(), v.size() * sizeof(T),
                    cudaMemcpyHostToDevice, stream);
    _size = v.size();
  }

  template <typename CPUAlloc = allocator<T>>
  void push_to(vector<T, CPUAlloc> &v, cudaStream_t stream = 0) const {
    v.resize(_size);
    cudaMemcpyAsync(v.data(), data_, _size * sizeof(T), cudaMemcpyDeviceToHost,
                    stream);
  }

  void clear() {
    _size = 0;
  }

  void reserve(int new_capacity) {
    if (new_capacity <= _capacity)
      return;

    auto new_data = alloc.allocate(new_capacity);
    if (_size > 0)
      cudaMemcpy(new_data, data_, _size * sizeof(T), cudaMemcpyDeviceToDevice);
    if (_capacity > 0)
      alloc.deallocate(data_, _capacity);

    data_ = new_data;
    _capacity = new_capacity;
  }

  T *data() {
    return data_;
  }

  T const *data() const {
    return data_;
  }

  __host__ __device__ int size() const {
    return _size;
  }

  operator view<T>() {
    return view<T>(data_, size());
  }

  operator const_view<T>() const {
    return const_view<T>(data_, size());
  }

private:
  T *data_;
  int _size, _capacity;
  Alloc alloc;

  void destroy() {
    if (_capacity > 0) {
      alloc.deallocate(data_, _capacity);
      data_ = nullptr;
      _size = _capacity = 0;
    }
  }
};
} // namespace nitro::def