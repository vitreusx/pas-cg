#pragma once
#include "byte.h"
#include <vcl/vectorclass.h>

namespace nitro::bit {
inline void load(Vec8ib &data, byte const *src) {
  data.load_bits(*src);
}

inline void store(Vec8ib const &data, byte *dst) {
  *dst = to_bits(data);
}

inline void load(Vec8qb &data, byte const *src) {
  data.load_bits(*src);
}

inline void store(Vec8qb const &data, byte *dst) {
  *dst = to_bits(data);
}

inline void load(Vec16ib &data, byte const *src) {
  data.load_bits(*(uint16_t const *)src);
}

inline void store(Vec16ib &data, byte *dst) {
  *(uint16_t *)dst = to_bits(data);
}
} // namespace nitro::bit