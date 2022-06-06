#pragma once
#include "vfunc.h"

#define _INST_CTOR_TLIST(...) VFUNC_B(_INST_CTOR_TLIST_, __VA_ARGS__)
#define _INST_CTOR_TLIST_1(x0) typename U0
#define _INST_CTOR_TLIST_2(x0, ...) _INST_CTOR_TLIST_1(__VA_ARGS__), typename U1
#define _INST_CTOR_TLIST_3(x0, ...) _INST_CTOR_TLIST_2(__VA_ARGS__), typename U2
#define _INST_CTOR_TLIST_4(x0, ...) _INST_CTOR_TLIST_3(__VA_ARGS__), typename U3
#define _INST_CTOR_TLIST_5(x0, ...) _INST_CTOR_TLIST_4(__VA_ARGS__), typename U4
#define _INST_CTOR_TLIST_6(x0, ...) _INST_CTOR_TLIST_5(__VA_ARGS__), typename U5
#define _INST_CTOR_TLIST_7(x0, ...) _INST_CTOR_TLIST_6(__VA_ARGS__), typename U6
#define _INST_CTOR_TLIST_8(x0, ...) _INST_CTOR_TLIST_7(__VA_ARGS__), typename U7