#pragma once
#include "const_ref.h"
#include <iterator>

namespace nitro::bit {
class const_iterator {
public:
  inline explicit const_iterator(byte const *ptr, int offset = 0) {
    this->ptr = ptr + byte_offset(offset);
    mask = byte_mask(offset);
  }

  inline explicit const_iterator(const_ref const &ref_)
      : ptr{ref_.p}, mask{ref_.mask} {}

  inline const_ref operator*() const {
    return const_ref(ptr, mask);
  }

  inline auto &operator++() {
    mask <<= 1;
    if (!mask) {
      ++ptr;
      mask = byte_mask(0);
    }
    return *this;
  }

  inline auto operator++(int) {
    auto sav = *this;
    ++*this;
    return sav;
  }

  inline auto &operator+=(int offset) {
    auto ptr_shift = byte_offset(offset),
         mask_shift = offset - (int)num_bits * ptr_shift;
    ptr += ptr_shift;
    while (mask_shift > 0)
      ++*this;
    return *this;
  }

  inline auto operator+(int offset) const {
    auto res = *this;
    res += offset;
    return res;
  }

  inline auto &operator--() {
    mask >>= 1;
    if (!mask) {
      --ptr;
      mask = byte_mask(num_bits - 1);
    }
    return *this;
  }

  inline auto operator--(int) {
    auto sav = *this;
    --*this;
    return sav;
  }

  inline auto &operator-=(int offset) {
    auto ptr_shift = byte_offset(offset),
         mask_shift = offset - (int)num_bits * ptr_shift;
    ptr -= ptr_shift;
    while (mask_shift > 0)
      --*this;
    return *this;
  }

  inline auto operator-(int offset) const {
    auto res = *this;
    res -= offset;
    return res;
  }

  inline std::ptrdiff_t operator-(const_iterator const &other) {
    return ptr - other.ptr;
  }

  inline bool operator<(const_iterator const &other) const {
    return ptr < other.ptr;
  }

  inline bool operator<=(const_iterator const &other) const {
    return ptr <= other.ptr;
  }

  inline bool operator>(const_iterator const &other) const {
    return ptr > other.ptr;
  }

  inline bool operator>=(const_iterator const &other) const {
    return ptr >= other.ptr;
  }

  inline bool operator==(const_iterator const &other) const {
    return ptr == other.ptr;
  }

  inline bool operator!=(const_iterator const &other) const {
    return ptr != other.ptr;
  }

private:
  byte const *ptr;
  byte mask;
};
} // namespace nitro::bit

namespace std {
template <>
struct iterator_traits<nitro::bit::const_iterator> {
  using value_type = bool;
  using difference_type = std::ptrdiff_t;
  using reference = nitro::bit::const_ref;
  using iterator_category = std::random_access_iterator_tag;
};
} // namespace std