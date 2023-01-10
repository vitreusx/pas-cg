#pragma once
#include "field.h"
#include "vfunc.h"

#define _INST_FIELDS(...) VFUNC(_INST_FIELDS_, __VA_ARGS__)

#define _INST_FIELDS_1(f0) _INST_FIELD_X(f0, NAME_OF(f0))
#define _INST_FIELD_X(f, x) _INST_FIELD_X2(f, x)
#define _INST_FIELD_X2(f, x)                                                   \
  __host__ __device__ TYPE_OF(f) & x() { return _##x; }                        \
                                                                               \
  __host__ __device__ TYPE_OF(f) const &x() const { return _##x; }             \
                                                                               \
  TYPE_OF(f) _##x;

#define _INST_FIELDS_2(f0, ...) _INST_FIELDS_1(f0) _INST_FIELDS_1(__VA_ARGS__)
#define _INST_FIELDS_3(f0, ...)                                                \
  _INST_FIELDS_1(f0)                                                           \
  _INST_FIELDS_2(__VA_ARGS__)
#define _INST_FIELDS_4(f0, ...)                                                \
  _INST_FIELDS_1(f0)                                                           \
  _INST_FIELDS_3(__VA_ARGS__)
#define _INST_FIELDS_5(f0, ...)                                                \
  _INST_FIELDS_1(f0)                                                           \
  _INST_FIELDS_4(__VA_ARGS__)
#define _INST_FIELDS_6(f0, ...)                                                \
  _INST_FIELDS_1(f0)                                                           \
  _INST_FIELDS_5(__VA_ARGS__)
#define _INST_FIELDS_7(f0, ...)                                                \
  _INST_FIELDS_1(f0)                                                           \
  _INST_FIELDS_6(__VA_ARGS__)
#define _INST_FIELDS_8(f0, ...)                                                \
  _INST_FIELDS_1(f0)                                                           \
  _INST_FIELDS_7(__VA_ARGS__)
