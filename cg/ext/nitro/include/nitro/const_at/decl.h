#pragma once
#include "../indexed.h"
namespace nitro {

template <typename T> struct const_at_expr_impl;

template <typename T>
using const_at_expr = typename const_at_expr_impl<T>::type;
} // namespace nitro