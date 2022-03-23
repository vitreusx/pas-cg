#pragma once
#include "config.h"

#define NARG(...) __NARG_I_(__VA_ARGS__, __RSEQ_N())
#define __NARG_I_(...) __ARG_N(__VA_ARGS__)
#define __ARG_N(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14,   \
                _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26,    \
                _27, _28, _29, _30, _31, _32, _33, _34, _35, _36, _37, _38,    \
                _39, _40, _41, _42, _43, _44, _45, _46, _47, _48, _49, _50,    \
                _51, _52, _53, _54, _55, _56, _57, _58, _59, _60, _61, _62,    \
                _63, N, ...)                                                   \
  N

#define __RSEQ_N()                                                             \
  63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45,  \
      44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27,  \
      26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9,   \
      8, 7, 6, 5, 4, 3, 2, 1, 0

#define __VFUNC(name, n) name##n
#define _VFUNC(name, n) __VFUNC(name, n)
#define VFUNC(func, ...) _VFUNC(func, NARG(__VA_ARGS__))(__VA_ARGS__)

#define EXPR_FIELD(name)                                                       \
  decltype(auto)  name() {                                         \
    return static_cast<E &>(*this).name();                                     \
  }                                                                            \
  decltype(auto)  name() const {                                   \
    return static_cast<E const &>(*this).name();                               \
  }

#define EXPR_FIELDS(...) VFUNC(EXPR_FIELDS_, __VA_ARGS__)

#define EXPR_FIELDS_1(x1) EXPR_FIELD(x1)
#define EXPR_FIELDS_2(x1, ...)                                                 \
  EXPR_FIELD(x1)                                                               \
  EXPR_FIELDS_1(__VA_ARGS__)
#define EXPR_FIELDS_3(x1, ...)                                                 \
  EXPR_FIELD(x1)                                                               \
  EXPR_FIELDS_2(__VA_ARGS__)
#define EXPR_FIELDS_4(x1, ...)                                                 \
  EXPR_FIELD(x1)                                                               \
  EXPR_FIELDS_3(__VA_ARGS__)
#define EXPR_FIELDS_5(x1, ...)                                                 \
  EXPR_FIELD(x1)                                                               \
  EXPR_FIELDS_4(__VA_ARGS__)
#define EXPR_FIELDS_6(x1, ...)                                                 \
  EXPR_FIELD(x1)                                                               \
  EXPR_FIELDS_5(__VA_ARGS__)
#define EXPR_FIELDS_7(x1, ...)                                                 \
  EXPR_FIELD(x1)                                                               \
  EXPR_FIELDS_6(__VA_ARGS__)
#define EXPR_FIELDS_8(x1, ...)                                                 \
  EXPR_FIELD(x1)                                                               \
  EXPR_FIELDS_7(__VA_ARGS__)

#define GET_BODY(...) GET_BODY_IX(0, __VA_ARGS__)
#define GET_BODY_IX(...) VFUNC(GET_BODY_IX_, __VA_ARGS__)

#define GET_BODY_IX_2(I0, x1)                                                  \
  if constexpr (I == I0)                                                       \
    return x1();
#define GET_BODY_IX_3(I0, x1, ...)                                             \
  if constexpr (I == I0)                                                       \
    return x1();                                                               \
  else                                                                         \
    GET_BODY_IX_2(I0 + 1, __VA_ARGS__)
#define GET_BODY_IX_4(I0, x1, ...)                                             \
  if constexpr (I == I0)                                                       \
    return x1();                                                               \
  else                                                                         \
    GET_BODY_IX_3(I0 + 1, __VA_ARGS__)
#define GET_BODY_IX_5(I0, x1, ...)                                             \
  if constexpr (I == I0)                                                       \
    return x1();                                                               \
  else                                                                         \
    GET_BODY_IX_4(I0 + 1, __VA_ARGS__)
#define GET_BODY_IX_6(I0, x1, ...)                                             \
  if constexpr (I == I0)                                                       \
    return x1();                                                               \
  else                                                                         \
    GET_BODY_IX_5(I0 + 1, __VA_ARGS__)
#define GET_BODY_IX_7(I0, x1, ...)                                             \
  if constexpr (I == I0)                                                       \
    return x1();                                                               \
  else                                                                         \
    GET_BODY_IX_6(I0 + 1, __VA_ARGS__)
#define GET_BODY_IX_8(I0, x1, ...)                                             \
  if constexpr (I == I0)                                                       \
    return x1();                                                               \
  else                                                                         \
    GET_BODY_IX_7(I0 + 1, __VA_ARGS__)
#define GET_BODY_IX_9(I0, x1, ...)                                             \
  if constexpr (I == I0)                                                       \
    return x1();                                                               \
  else                                                                         \
    GET_BODY_IX_8(I0 + 1, __VA_ARGS__)

#define EXPR_GETS(...)                                                         \
  template <size_t I>  decltype(auto) get() {                      \
    GET_BODY(__VA_ARGS__)                                                      \
  }                                                                            \
  template <size_t I>  decltype(auto) get() const {                \
    GET_BODY(__VA_ARGS__)                                                      \
  }

#define EXPR_BODY(...)                                                         \
  EXPR_FIELDS(__VA_ARGS__)                                                     \
  EXPR_GETS(__VA_ARGS__)

#define AUTO_EXPR_FIELD(I, name)                                               \
   decltype(auto) name() {                                         \
    return static_cast<E &>(*this).template get<I>();                          \
  }                                                                            \
   decltype(auto) name() const {                                   \
    return static_cast<E const &>(*this).template get<I>();                    \
  }

#define AUTO_EXPR_FIELDS(...) AUTO_EXPR_FIELDS_IX(0, __VA_ARGS__)
#define AUTO_EXPR_FIELDS_IX(...) VFUNC(AUTO_EXPR_FIELDS_IX_, __VA_ARGS__)

#define AUTO_EXPR_FIELDS_IX_2(I, x1) AUTO_EXPR_FIELD(I, x1)
#define AUTO_EXPR_FIELDS_IX_3(I, x1, ...)                                      \
  AUTO_EXPR_FIELD(I, x1)                                                       \
  AUTO_EXPR_FIELDS_IX_2(I + 1, __VA_ARGS__)
#define AUTO_EXPR_FIELDS_IX_4(I, x1, ...)                                      \
  AUTO_EXPR_FIELD(I, x1)                                                       \
  AUTO_EXPR_FIELDS_IX_3(I + 1, __VA_ARGS__)
#define AUTO_EXPR_FIELDS_IX_5(I, x1, ...)                                      \
  AUTO_EXPR_FIELD(I, x1)                                                       \
  AUTO_EXPR_FIELDS_IX_4(I + 1, __VA_ARGS__)
#define AUTO_EXPR_FIELDS_IX_6(I, x1, ...)                                      \
  AUTO_EXPR_FIELD(I, x1)                                                       \
  AUTO_EXPR_FIELDS_IX_5(I + 1, __VA_ARGS__)
#define AUTO_EXPR_FIELDS_IX_7(I, x1, ...)                                      \
  AUTO_EXPR_FIELD(I, x1)                                                       \
  AUTO_EXPR_FIELDS_IX_6(I + 1, __VA_ARGS__)
#define AUTO_EXPR_FIELDS_IX_8(I, x1, ...)                                      \
  AUTO_EXPR_FIELD(I, x1)                                                       \
  AUTO_EXPR_FIELDS_IX_7(I + 1, __VA_ARGS__)
#define AUTO_EXPR_FIELDS_IX_9(I, x1, ...)                                      \
  AUTO_EXPR_FIELD(I, x1)                                                       \
  AUTO_EXPR_FIELDS_IX_8(I + 1, __VA_ARGS__)

#define AUTO_EXPR_GETS()                                                       \
  template <size_t I>  decltype(auto) get() {                      \
    return static_cast<E &>(*this).template get<I>();                          \
  }                                                                            \
  template <size_t I>  decltype(auto) get() const {                \
    return static_cast<E const &>(*this).template get<I>();                    \
  }

#define AUTO_EXPR_BODY(...)                                                    \
  AUTO_EXPR_FIELDS(__VA_ARGS__)                                                \
  AUTO_EXPR_GETS()
