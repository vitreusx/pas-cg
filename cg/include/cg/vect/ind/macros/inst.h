#pragma once
#include "../type_list.h"
#include "inst_ctor.h"
#include "inst_expr_ctor.h"
#include "inst_fields.h"
#include "inst_subtypes.h"

#define INST(cls, ...)                                                         \
  _INST_CTOR(cls, __VA_ARGS__)                                                 \
  _INST_EXPR_CTOR(cls, __VA_ARGS__)                                            \
  _INST_FIELDS(__VA_ARGS__)                                                    \
                                                                               \
  __host__ __device__ static auto Subtypes() {                                 \
    return nitro::ind::type_list<_INST_SUBTYPES(__VA_ARGS__)>{};               \
  }                                                                            \
                                                                               \
  template <typename E>                                                        \
  __host__ __device__ static auto Expr()->cls##_expr<E>;
