#pragma once
#include "vfunc.h"

#define _IEXPR_NAMED(...) VFUNC(_IEXPR_NAMED_, 0, __VA_ARGS__)

#define _IEXPR_NAMED_2(I_x, x)                                                 \
  __host__ __device__ decltype(auto) x() {                                     \
    return (static_cast<E &>(*this).template get<I_x>());                      \
  }                                                                            \
                                                                               \
  __host__ __device__ decltype(auto) x() const {                               \
    return (static_cast<E const &>(*this).template get<I_x>());                \
  }

#define _IEXPR_NAMED_3(I_0, x_0, ...)                                          \
  _IEXPR_NAMED_2(I_0, x_0)                                                     \
  _IEXPR_NAMED_2(I_0 + 1, __VA_ARGS__)

#define _IEXPR_NAMED_4(I_0, x_0, ...)                                          \
  _IEXPR_NAMED_2(I_0, x_0)                                                     \
  _IEXPR_NAMED_3(I_0 + 1, __VA_ARGS__)

#define _IEXPR_NAMED_5(I_0, x_0, ...)                                          \
  _IEXPR_NAMED_2(I_0, x_0)                                                     \
  _IEXPR_NAMED_4(I_0 + 1, __VA_ARGS__)

#define _IEXPR_NAMED_6(I_0, x_0, ...)                                          \
  _IEXPR_NAMED_2(I_0, x_0)                                                     \
  _IEXPR_NAMED_5(I_0 + 1, __VA_ARGS__)

#define _IEXPR_NAMED_7(I_0, x_0, ...)                                          \
  _IEXPR_NAMED_2(I_0, x_0)                                                     \
  _IEXPR_NAMED_6(I_0 + 1, __VA_ARGS__)

#define _IEXPR_NAMED_8(I_0, x_0, ...)                                          \
  _IEXPR_NAMED_2(I_0, x_0)                                                     \
  _IEXPR_NAMED_7(I_0 + 1, __VA_ARGS__)

#define _IEXPR_NAMED_9(I_0, x_0, ...)                                          \
  _IEXPR_NAMED_2(I_0, x_0)                                                     \
  _IEXPR_NAMED_8(I_0 + 1, __VA_ARGS__)
