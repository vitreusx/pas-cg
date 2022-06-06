#pragma once
#include "narg.h"

#define __VFUNC(name, n) name##n
#define _VFUNC(name, n) __VFUNC(name, n)
#define VFUNC(func, ...) _VFUNC(func, NARG(__VA_ARGS__))(__VA_ARGS__)

#define __VFUNC_B(name, n) name##n
#define _VFUNC_B(name, n) __VFUNC_B(name, n)
#define VFUNC_B(func, ...) _VFUNC_B(func, NARG(__VA_ARGS__))(__VA_ARGS__)

#define __VFUNC_C(name, n) name##n
#define _VFUNC_C(name, n) __VFUNC_C(name, n)
#define VFUNC_C(func, ...) _VFUNC_C(func, NARG(__VA_ARGS__))(__VA_ARGS__)
