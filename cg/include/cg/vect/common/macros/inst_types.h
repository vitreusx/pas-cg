#pragma once
#include "field.h"
#include "vfunc.h"

#define _INST_TYPES(...) VFUNC(_INST_TYPES_, __VA_ARGS__)
#define _INST_TYPES_1(f0) TYPE_OF(f0)
#define _INST_TYPES_2(f0, ...) _INST_TYPES_1(f0), _INST_TYPES_1(__VA_ARGS__)
#define _INST_TYPES_3(f0, ...) _INST_TYPES_1(f0), _INST_TYPES_2(__VA_ARGS__)
#define _INST_TYPES_4(f0, ...) _INST_TYPES_1(f0), _INST_TYPES_3(__VA_ARGS__)
#define _INST_TYPES_5(f0, ...) _INST_TYPES_1(f0), _INST_TYPES_4(__VA_ARGS__)
#define _INST_TYPES_6(f0, ...) _INST_TYPES_1(f0), _INST_TYPES_5(__VA_ARGS__)
#define _INST_TYPES_7(f0, ...) _INST_TYPES_1(f0), _INST_TYPES_6(__VA_ARGS__)
#define _INST_TYPES_8(f0, ...) _INST_TYPES_1(f0), _INST_TYPES_7(__VA_ARGS__)
