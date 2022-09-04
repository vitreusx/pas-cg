#pragma once
#include "vfunc.h"

#define _INST_CTOR_SLIST(...) VFUNC_B(_INST_CTOR_SLIST_, __VA_ARGS__)
#define _INST_CTOR_SLIST_1(x0) _##x0(x0)
#define _INST_CTOR_SLIST_2(x0, ...)                                            \
  _INST_CTOR_SLIST_1(x0), _INST_CTOR_SLIST_1(__VA_ARGS__)
#define _INST_CTOR_SLIST_3(x0, ...)                                            \
  _INST_CTOR_SLIST_1(x0), _INST_CTOR_SLIST_2(__VA_ARGS__)
#define _INST_CTOR_SLIST_4(x0, ...)                                            \
  _INST_CTOR_SLIST_1(x0), _INST_CTOR_SLIST_3(__VA_ARGS__)
#define _INST_CTOR_SLIST_5(x0, ...)                                            \
  _INST_CTOR_SLIST_1(x0), _INST_CTOR_SLIST_4(__VA_ARGS__)
#define _INST_CTOR_SLIST_6(x0, ...)                                            \
  _INST_CTOR_SLIST_1(x0), _INST_CTOR_SLIST_5(__VA_ARGS__)
#define _INST_CTOR_SLIST_7(x0, ...)                                            \
  _INST_CTOR_SLIST_1(x0), _INST_CTOR_SLIST_6(__VA_ARGS__)
#define _INST_CTOR_SLIST_8(x0, ...)                                            \
  _INST_CTOR_SLIST_1(x0), _INST_CTOR_SLIST_7(__VA_ARGS__)
