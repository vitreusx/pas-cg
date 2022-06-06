#pragma once

#define _EXPR_GET_BODY(...) VFUNC(_EXPR_GET_BODY_, 0, __VA_ARGS__)

#define _EXPR_GET_BODY_3(I_x, qual, x)                                         \
  if constexpr (I == I_x)                                                      \
    return (static_cast<E qual &>(*this).x());

#define _EXPR_GET_BODY_4(I_0, qual, x0, ...)                                   \
  _EXPR_GET_BODY_3(I_0, qual, x0)                                              \
  else _EXPR_GET_BODY_3(I_0 + 1, qual, __VA_ARGS__)

#define _EXPR_GET_BODY_5(I_0, qual, x0, ...)                                   \
  _EXPR_GET_BODY_3(I_0, qual, x0)                                              \
  else _EXPR_GET_BODY_4(I_0 + 1, qual, __VA_ARGS__)

#define _EXPR_GET_BODY_6(I_0, qual, x0, ...)                                   \
  _EXPR_GET_BODY_3(I_0, qual, x0)                                              \
  else _EXPR_GET_BODY_5(I_0 + 1, qual, __VA_ARGS__)

#define _EXPR_GET_BODY_7(I_0, qual, x0, ...)                                   \
  _EXPR_GET_BODY_3(I_0, qual, x0)                                              \
  else _EXPR_GET_BODY_6(I_0 + 1, qual, __VA_ARGS__)

#define _EXPR_GET_BODY_8(I_0, qual, x0, ...)                                   \
  _EXPR_GET_BODY_3(I_0, qual, x0)                                              \
  else _EXPR_GET_BODY_7(I_0 + 1, qual, __VA_ARGS__)

#define _EXPR_GET_BODY_9(I_0, qual, x0, ...)                                   \
  _EXPR_GET_BODY_3(I_0, qual, x0)                                              \
  else _EXPR_GET_BODY_8(I_0 + 1, qual, __VA_ARGS__)

#define _EXPR_GET_BODY_10(I_0, qual, x0, ...)                                  \
  _EXPR_GET_BODY_3(I_0, qual, x0)                                              \
  else _EXPR_GET_BODY_9(I_0 + 1, qual, __VA_ARGS__)

#define _EXPR_GET_BODY_11(I_0, qual, x0, ...)                                  \
  _EXPR_GET_BODY_3(I_0, qual, x0)                                              \
  else _EXPR_GET_BODY_10(I_0 + 1, qual, __VA_ARGS__)
