#pragma once
#include "inst_expr_ctor_slist.h"

#define _INST_EXPR_CTOR(cls, ...)                                              \
  template <typename E>                                                        \
  __host__ __device__ cls(cls##_expr<E> const &e)                              \
      : _INST_EXPR_CTOR_SLIST(__VA_ARGS__) {}
