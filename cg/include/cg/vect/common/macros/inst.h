#pragma once
#include "../type_list.h"
#include "inst_ctor.h"
#include "inst_expr_ctor.h"
#include "inst_fields.h"
#include "inst_types.h"

#define INST(cls, ...)                                                         \
  _INST_CTOR(cls, __VA_ARGS__)                                                 \
  _INST_EXPR_CTOR(cls, __VA_ARGS__)                                            \
  _INST_FIELDS(__VA_ARGS__)                                                    \
                                                                               \
  static auto Types() {                                                        \
    return type_list<_INST_TYPES(__VA_ARGS__)>{};                              \
  }
