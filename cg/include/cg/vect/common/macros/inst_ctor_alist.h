#pragma once
#include "vfunc.h"

#define _INST_CTOR_ALIST(...) VFUNC_B(_INST_CTOR_ALIST_, __VA_ARGS__)
#define _INST_CTOR_ALIST_1(x0) U0 &&x0
#define _INST_CTOR_ALIST_2(x0, x1) U0 &&x0, U1 &&x1
#define _INST_CTOR_ALIST_3(x0, x1, x2) U0 &&x0, U1 &&x1, U2 &&x2
#define _INST_CTOR_ALIST_4(x0, x1, x2, x3) U0 &&x0, U1 &&x1, U2 &&x2, U3 &&x3
#define _INST_CTOR_ALIST_5(x0, x1, x2, x3, x4)                                 \
  U0 &&x0, U1 &&x1, U2 &&x2, U3 &&x3, U4 &&x4
#define _INST_CTOR_ALIST_6(x0, x1, x2, x3, x4, x5)                             \
  U0 &&x0, U1 &&x1, U2 &&x2, U3 &&x3, U4 &&x4, U5 &&x5
#define _INST_CTOR_ALIST_7(x0, x1, x2, x3, x4, x5, x6)                         \
  U0 &&x0, U1 &&x1, U2 &&x2, U3 &&x3, U4 &&x4, U5 &&x5, U6 &&x6
#define _INST_CTOR_ALIST_8(x0, x1, x2, x3, x4, x5, x6, x7)                     \
  U0 &&x0, U1 &&x1, U2 &&x2, U3 &&x3, U4 &&x4, U5 &&x5, U6 &&x6, U7 &&x7
