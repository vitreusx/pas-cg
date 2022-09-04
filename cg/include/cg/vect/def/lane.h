#pragma once
#include <cstdint>
#include <type_traits>
#include <vcl/vectorclass.h>

namespace nitro::def {
template <std::size_t Nb>
struct _uint_from_size_impl;

template <std::size_t Nb>
using uint_from_size_t = typename _uint_from_size_impl<Nb>::type;

#define IMPL(Nb, T)                                                            \
  template <>                                                                  \
  struct _uint_from_size_impl<Nb> {                                            \
    using type = T;                                                            \
  }

IMPL(8, uint8_t);
IMPL(16, uint16_t);
IMPL(32, uint32_t);
IMPL(64, uint64_t);

#undef IMPL

template <typename T>
using uint_repr = uint_from_size_t<8 * sizeof(T)>;

template <typename T, std::size_t N, std::size_t W>
struct _lane_impl {
  using type = typename _lane_impl<uint_repr<T>, N, W>::type;
};

template <typename T, std::size_t N, std::size_t W>
using lane = typename _lane_impl<T, N, W>::type;

#define IMPL(T, N, W, L)                                                       \
  template <>                                                                  \
  struct _lane_impl<T, N, W> {                                                 \
    using type = L;                                                            \
  }

IMPL(int8_t, 8, 256, Vec8i);
IMPL(uint8_t, 8, 256, Vec8ui);
IMPL(int16_t, 8, 256, Vec8i);
IMPL(uint16_t, 8, 256, Vec8ui);
IMPL(int32_t, 8, 256, Vec8i);
IMPL(uint32_t, 8, 256, Vec8ui);
IMPL(float, 8, 256, Vec8f);
IMPL(bool, 8, 256, Vec8ib);

IMPL(int8_t, 8, 512, Vec8q);
IMPL(uint8_t, 8, 512, Vec8uq);
IMPL(int16_t, 8, 512, Vec8q);
IMPL(uint16_t, 8, 512, Vec8uq);
IMPL(int32_t, 8, 512, Vec8q);
IMPL(uint32_t, 8, 512, Vec8uq);
IMPL(int64_t, 8, 512, Vec8q);
IMPL(uint64_t, 8, 512, Vec8uq);
IMPL(float, 8, 512, Vec8d);
IMPL(double, 8, 512, Vec8d);
IMPL(bool, 8, 512, Vec8qb);

IMPL(int8_t, 16, 512, Vec16i);
IMPL(uint8_t, 16, 512, Vec16ui);
IMPL(int16_t, 16, 512, Vec16i);
IMPL(uint16_t, 16, 512, Vec16ui);
IMPL(int32_t, 16, 512, Vec16i);
IMPL(uint32_t, 16, 512, Vec16ui);
IMPL(float, 16, 512, Vec16f);
IMPL(bool, 16, 512, Vec16ib);

#undef IMPL

} // namespace nitro