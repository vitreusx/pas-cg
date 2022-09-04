#pragma once
#include "byte.h"

namespace nitro::bit {
struct ref {
public:
  inline explicit ref(byte *p, byte mask) : p{p}, mask{mask} {}

  friend void swap(ref const &x, ref const &y) {
    bool tmp = y;
    y = x;
    x = tmp;
  }

  inline operator bool() const {
    return (*p) & mask;
  }

  inline bool operator!() const {
    return !(bool)*this;
  }

  template <typename U>
  ref const &operator=(U const &val) const {
    if (val)
      *p |= mask;
    else
      *p &= (~mask);
    return *this;
  }

  byte *p, mask;
};
} // namespace nitro::bit