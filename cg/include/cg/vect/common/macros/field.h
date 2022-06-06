#pragma once
#include "vfunc.h"

#define FIELD(...) FIELD(__VA_ARGS__)

#define TYPE_OF(f) TYPE_OF_##f
#define TYPE_OF_FIELD(...) VFUNC_C(TYPE_OF_FIELD_, __VA_ARGS__)
#define TYPE_OF_FIELD_2(T, x) T
#define TYPE_OF_FIELD_3(T1, ...) T1, TYPE_OF_FIELD_2(__VA_ARGS__)
#define TYPE_OF_FIELD_4(T1, ...) T1, TYPE_OF_FIELD_3(__VA_ARGS__)
#define TYPE_OF_FIELD_5(T1, ...) T1, TYPE_OF_FIELD_4(__VA_ARGS__)
#define TYPE_OF_FIELD_6(T1, ...) T1, TYPE_OF_FIELD_5(__VA_ARGS__)
#define TYPE_OF_FIELD_7(T1, ...) T1, TYPE_OF_FIELD_6(__VA_ARGS__)
#define TYPE_OF_FIELD_8(T1, ...) T1, TYPE_OF_FIELD_7(__VA_ARGS__)
#define TYPE_OF_FIELD_9(T1, ...) T1, TYPE_OF_FIELD_8(__VA_ARGS__)

#define NAME_OF(f) NAME_OF_##f
#define NAME_OF_FIELD(...) VFUNC_C(NAME_OF_FIELD_, __VA_ARGS__)
#define NAME_OF_FIELD_2(T, x) x
#define NAME_OF_FIELD_3(T, ...) NAME_OF_FIELD_2(__VA_ARGS__)
#define NAME_OF_FIELD_4(T, ...) NAME_OF_FIELD_3(__VA_ARGS__)
#define NAME_OF_FIELD_5(T, ...) NAME_OF_FIELD_4(__VA_ARGS__)
#define NAME_OF_FIELD_6(T, ...) NAME_OF_FIELD_5(__VA_ARGS__)
#define NAME_OF_FIELD_7(T, ...) NAME_OF_FIELD_6(__VA_ARGS__)
#define NAME_OF_FIELD_8(T, ...) NAME_OF_FIELD_7(__VA_ARGS__)
#define NAME_OF_FIELD_9(T, ...) NAME_OF_FIELD_8(__VA_ARGS__)