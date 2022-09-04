#pragma once
#include <cstdint>

namespace nitro::bit {
using byte = uint8_t;
inline constexpr std::size_t num_bits = 8 * sizeof(byte);

template <typename Int>
inline byte byte_mask(Int offset) {
  return (byte)1 << (offset % num_bits);
}

template <typename Int>
inline Int byte_offset(Int offset) {
  return offset / num_bits;
}

template <typename Int>
inline Int req_bytes(Int n) {
  return (n + num_bits - 1) / num_bits;
}
} // namespace nitro::bit