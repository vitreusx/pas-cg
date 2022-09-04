#pragma once
#include "field.h"
#include "vfunc.h"

#define _INST_SUBTYPES(...) VFUNC(_INST_SUBTYPES_, __VA_ARGS__)
#define _INST_SUBTYPES_1(f0) TYPE_OF(f0)
#define _INST_SUBTYPES_2(f0, ...)                                              \
  _INST_SUBTYPES_1(f0), _INST_SUBTYPES_1(__VA_ARGS__)
#define _INST_SUBTYPES_3(f0, ...)                                              \
  _INST_SUBTYPES_1(f0), _INST_SUBTYPES_2(__VA_ARGS__)
#define _INST_SUBTYPES_4(f0, ...)                                              \
  _INST_SUBTYPES_1(f0), _INST_SUBTYPES_3(__VA_ARGS__)
#define _INST_SUBTYPES_5(f0, ...)                                              \
  _INST_SUBTYPES_1(f0), _INST_SUBTYPES_4(__VA_ARGS__)
#define _INST_SUBTYPES_6(f0, ...)                                              \
  _INST_SUBTYPES_1(f0), _INST_SUBTYPES_5(__VA_ARGS__)
#define _INST_SUBTYPES_7(f0, ...)                                              \
  _INST_SUBTYPES_1(f0), _INST_SUBTYPES_6(__VA_ARGS__)
#define _INST_SUBTYPES_8(f0, ...)                                              \
  _INST_SUBTYPES_1(f0), _INST_SUBTYPES_7(__VA_ARGS__)
