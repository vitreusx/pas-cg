#pragma once
#include "vfunc.h"

#define _INST_EXPR_CTOR_SLIST(...) VFUNC(_INST_EXPR_CTOR_SLIST_, __VA_ARGS__)

#define _INST_EXPR_CTOR_SLIST_1(f0) _INST_EXPR_CTOR_SLIST_X(NAME_OF(f0))
#define _INST_EXPR_CTOR_SLIST_X(x0) _INST_EXPR_CTOR_SLIST_Y(x0)
#define _INST_EXPR_CTOR_SLIST_Y(x0) _##x0(e.x0())

#define _INST_EXPR_CTOR_SLIST_2(f0, ...)                                       \
  _INST_EXPR_CTOR_SLIST_1(f0), _INST_EXPR_CTOR_SLIST_1(__VA_ARGS__)
#define _INST_EXPR_CTOR_SLIST_3(f0, ...)                                       \
  _INST_EXPR_CTOR_SLIST_1(f0), _INST_EXPR_CTOR_SLIST_2(__VA_ARGS__)
#define _INST_EXPR_CTOR_SLIST_4(f0, ...)                                       \
  _INST_EXPR_CTOR_SLIST_1(f0), _INST_EXPR_CTOR_SLIST_3(__VA_ARGS__)
#define _INST_EXPR_CTOR_SLIST_5(f0, ...)                                       \
  _INST_EXPR_CTOR_SLIST_1(f0), _INST_EXPR_CTOR_SLIST_4(__VA_ARGS__)
#define _INST_EXPR_CTOR_SLIST_6(f0, ...)                                       \
  _INST_EXPR_CTOR_SLIST_1(f0), _INST_EXPR_CTOR_SLIST_5(__VA_ARGS__)
#define _INST_EXPR_CTOR_SLIST_7(f0, ...)                                       \
  _INST_EXPR_CTOR_SLIST_1(f0), _INST_EXPR_CTOR_SLIST_6(__VA_ARGS__)
#define _INST_EXPR_CTOR_SLIST_8(f0, ...)                                       \
  _INST_EXPR_CTOR_SLIST_1(f0), _INST_EXPR_CTOR_SLIST_7(__VA_ARGS__)
