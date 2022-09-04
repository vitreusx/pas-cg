#pragma once
#include "inst_ctor_alist.h"
#include "inst_ctor_slist.h"
#include "inst_ctor_tlist.h"
#include "vfunc.h"

#define _INST_CTOR(cls, ...) VFUNC(_INST_CTOR_, cls, __VA_ARGS__)
#define _INST_CTOR_2(cls, f0) _INST_CTOR_REAL(cls, NAME_OF(f0))
#define _INST_CTOR_3(cls, f0, f1) _INST_CTOR_REAL(cls, NAME_OF(f0), NAME_OF(f1))
#define _INST_CTOR_4(cls, f0, f1, f2)                                          \
  _INST_CTOR_REAL(cls, NAME_OF(f0), NAME_OF(f1), NAME_OF(f2))
#define _INST_CTOR_5(cls, f0, f1, f2, f3)                                      \
  _INST_CTOR_REAL(cls, NAME_OF(f0), NAME_OF(f1), NAME_OF(f2), NAME_OF(f3))
#define _INST_CTOR_6(cls, f0, f1, f2, f3, f4)                                  \
  _INST_CTOR_REAL(cls, NAME_OF(f0), NAME_OF(f1), NAME_OF(f2), NAME_OF(f3),     \
                  NAME_OF(f4))
#define _INST_CTOR_7(cls, f0, f1, f2, f3, f4, f5)                              \
  _INST_CTOR_REAL(cls, NAME_OF(f0), NAME_OF(f1), NAME_OF(f2), NAME_OF(f3),     \
                  NAME_OF(f4), NAME_OF(f5))
#define _INST_CTOR_8(cls, f0, f1, f2, f3, f4, f5, f6)                          \
  _INST_CTOR_REAL(cls, NAME_OF(f0), NAME_OF(f1), NAME_OF(f2), NAME_OF(f3),     \
                  NAME_OF(f4), NAME_OF(f5), NAME_OF(f6))
#define _INST_CTOR_9(cls, f0, f1, f2, f3, f4, f5, f6, f7)                      \
  _INST_CTOR_REAL(cls, NAME_OF(f0), NAME_OF(f1), NAME_OF(f2), NAME_OF(f3),     \
                  NAME_OF(f4), NAME_OF(f5), NAME_OF(f6), NAME_OF(f7))

#define _INST_CTOR_REAL(cls, ...)                                              \
  template <_INST_CTOR_TLIST(__VA_ARGS__)>                                     \
  cls(_INST_CTOR_ALIST(__VA_ARGS__)) : _INST_CTOR_SLIST(__VA_ARGS__) {}
