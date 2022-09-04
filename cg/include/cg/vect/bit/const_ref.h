#pragma once
#include "byte.h"

namespace nitro::bit {
struct const_ref {
public:
  inline explicit const_ref(byte const *p, byte mask) : p{p}, mask{mask} {}

  template <typename U>
  const_ref &operator=(U const &val) const = delete;

  inline operator bool() const {
    return (*p) & mask;
  }

  inline bool operator!() const {
    return !(bool)*this;
  }

  byte const *p;
  byte mask;
};
} // namespace nitro::bit