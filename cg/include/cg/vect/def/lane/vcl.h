#pragma once
#include "type_traits.h"
#include <tuple>
#include <type_traits>
#include <vcl/vectorclass.h>

namespace nitro::def {
template <typename T, std::size_t N, std::size_t W>
struct _vcl_lane;

template <typename T, std::size_t N, std::size_t W>
using vcl_lane = typename _vcl_lane<T, N, W>::type;

#define IMPL(T, N, W, L)                                                       \
  template <>                                                                  \
  struct _vcl_lane<T, N, W> {                                                  \
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

IMPL(int8_t, 4, 256, Vec4q);
IMPL(uint8_t, 4, 256, Vec4uq);
IMPL(int16_t, 4, 256, Vec4q);
IMPL(uint16_t, 4, 256, Vec4uq);
IMPL(int32_t, 4, 256, Vec4q);
IMPL(uint32_t, 4, 256, Vec4uq);
IMPL(int64_t, 4, 256, Vec4q);
IMPL(uint64_t, 4, 256, Vec4uq);
IMPL(float, 4, 256, Vec4d);
IMPL(double, 4, 256, Vec4d);
IMPL(bool, 4, 256, Vec4qb);

IMPL(int8_t, 16, 512, Vec16i);
IMPL(uint8_t, 16, 512, Vec16ui);
IMPL(int16_t, 16, 512, Vec16i);
IMPL(uint16_t, 16, 512, Vec16ui);
IMPL(int32_t, 16, 512, Vec16i);
IMPL(uint32_t, 16, 512, Vec16ui);
IMPL(float, 16, 512, Vec16f);
IMPL(bool, 16, 512, Vec16ib);

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

#undef IMPL

template <typename T>
struct is_vcl_lane : std::false_type {};

template <typename T>
inline constexpr bool is_vcl_lane_v = is_vcl_lane<T>::value;

template <typename T>
struct vcl_lane_type {
  using type = std::tuple_element_t<
      T::elementtype(),
      std::tuple<void *, bool, bool, bool, int8_t, uint8_t, int16_t, uint16_t,
                 int32_t, uint32_t, int64_t, uint64_t, void *, void *, void *,
                 void *, float, double>>;
};

#define VCL_LANE(T)                                                            \
  template <>                                                                  \
  struct is_lane_like<T> : std::true_type {};                                  \
                                                                               \
  template <>                                                                  \
  struct is_vcl_lane<T> : std::true_type {};                                   \
                                                                               \
  template <>                                                                  \
  struct lane_type<T> {                                                        \
    using type = typename vcl_lane_type<T>::type;                              \
  };                                                                           \
                                                                               \
  template <>                                                                  \
  struct lane_size<T> {                                                        \
    static constexpr std::size_t value = T::size();                            \
  }

VCL_LANE(Vec8i);
VCL_LANE(Vec8ui);
VCL_LANE(Vec8f);
VCL_LANE(Vec8ib);

VCL_LANE(Vec4q);
VCL_LANE(Vec4uq);
VCL_LANE(Vec4d);
VCL_LANE(Vec4qb);

VCL_LANE(Vec8q);
VCL_LANE(Vec8uq);
VCL_LANE(Vec8d);
#if INSTRSET < 10
VCL_LANE(Vec8qb);
#endif

VCL_LANE(Vec16i);
VCL_LANE(Vec16ui);
VCL_LANE(Vec16f);
VCL_LANE(Vec16ib);

#undef VCL_LANE

} // namespace nitro::def