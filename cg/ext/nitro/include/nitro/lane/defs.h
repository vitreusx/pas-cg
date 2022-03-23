#pragma once
#include <cstddef>
#include <cstdint>
#include <vcl/vectorclass.h>

namespace nitro {
template <typename T, size_t N> struct lane_impl;
template <typename T, size_t N> using lane = typename lane_impl<T, N>::type;

#define LANE_IMPL(T, N, V)                                                     \
  template <> struct lane_impl<T, N> { using type = V; }

LANE_IMPL(int8_t, 16, Vec16c);
LANE_IMPL(uint8_t, 16, Vec16uc);
LANE_IMPL(int16_t, 8, Vec8s);
LANE_IMPL(uint16_t, 8, Vec8us);
LANE_IMPL(int32_t, 4, Vec4i);
LANE_IMPL(uint32_t, 4, Vec4ui);
LANE_IMPL(int64_t, 2, Vec2q);
LANE_IMPL(uint64_t, 2, Vec2uq);
LANE_IMPL(float, 4, Vec4f);
LANE_IMPL(double, 2, Vec2d);

LANE_IMPL(int8_t, 32, Vec32c);
LANE_IMPL(uint8_t, 32, Vec32uc);
LANE_IMPL(int16_t, 16, Vec16s);
LANE_IMPL(uint16_t, 16, Vec16us);
LANE_IMPL(int32_t, 8, Vec8i);
LANE_IMPL(uint32_t, 8, Vec8ui);
LANE_IMPL(int64_t, 4, Vec4q);
LANE_IMPL(uint64_t, 4, Vec4uq);
LANE_IMPL(float, 8, Vec8f);
LANE_IMPL(double, 4, Vec4d);

LANE_IMPL(int8_t, 64, Vec64c);
LANE_IMPL(uint8_t, 64, Vec64uc);
LANE_IMPL(int16_t, 32, Vec32s);
LANE_IMPL(uint16_t, 32, Vec32us);
LANE_IMPL(int32_t, 16, Vec16i);
LANE_IMPL(uint32_t, 16, Vec16ui);
LANE_IMPL(int64_t, 8, Vec8q);
LANE_IMPL(uint64_t, 8, Vec8uq);
LANE_IMPL(float, 16, Vec16f);
LANE_IMPL(double, 8, Vec8d);

#undef LANE_IMPL

template <typename T> struct mask_impl;
template <typename T> using mask = typename mask_impl<T>::type;

#define MASK_IMPL(T, M)                                                        \
  template <> struct mask_impl<T> { using type = M; }

#define MASK2_IMPL(T1, T2, M)                                                  \
  MASK_IMPL(T1, M);                                                            \
  MASK_IMPL(T2, M)

MASK2_IMPL(Vec16c, Vec16uc, Vec16cb);
MASK2_IMPL(Vec8s, Vec8us, Vec8sb);
MASK2_IMPL(Vec4i, Vec4ui, Vec4ib);
MASK2_IMPL(Vec2q, Vec2uq, Vec2qb);
MASK_IMPL(Vec4f, Vec4fb);
MASK_IMPL(Vec2d, Vec2db);
MASK_IMPL(Vec8f, Vec8fb);
MASK_IMPL(Vec4d, Vec4db);

MASK2_IMPL(Vec32c, Vec32uc, Vec32cb);
MASK2_IMPL(Vec16s, Vec16us, Vec16sb);
MASK2_IMPL(Vec8i, Vec8ui, Vec8ib);
MASK2_IMPL(Vec4q, Vec4uq, Vec4qb);

MASK2_IMPL(Vec64c, Vec64uc, Vec64cb);
MASK2_IMPL(Vec32s, Vec32us, Vec32sb);
MASK2_IMPL(Vec16i, Vec16ui, Vec16ib);
MASK2_IMPL(Vec8q, Vec8uq, Vec8qb);
MASK_IMPL(Vec16f, Vec16fb);
MASK_IMPL(Vec8d, Vec8db);

#undef MASK2_IMPL
#undef MASK_IMPl
}
