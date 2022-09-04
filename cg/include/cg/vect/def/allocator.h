#pragma once
#include <cstdint>
#include <cstdlib>
#include <memory>

namespace nitro::def {
template <typename T>
struct allocator {
  T *allocate(std::size_t n, const void * = nullptr) {
    std::size_t alignment = 512 >> 3, nbytes = n * sizeof(T);
    nbytes = alignment * ((nbytes + alignment - 1) / alignment);
    return (T *)std::aligned_alloc(alignment, nbytes);
    //    return (T *)std::malloc(n * sizeof(T));
    //    return alloc.allocate(n);
  }

  void deallocate(T *p, std::size_t) {
    free((void *)p);
    //    alloc.deallocate(p, n);
  }

  //  std::allocator<T> alloc;
};
} // namespace nitro::def