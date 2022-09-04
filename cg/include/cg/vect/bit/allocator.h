#pragma once
#include "../def/allocator.h"
#include "byte.h"

namespace nitro::bit {
struct allocator {
  inline byte *allocate(std::size_t n, const void *hint = nullptr) {
    return alloc.allocate(req_bytes(n), hint);
  }

  inline void deallocate(byte *p, std::size_t n) {
    alloc.deallocate(p, req_bytes(n));
  }

  def::allocator<byte> alloc;
};
} // namespace nitro::bit