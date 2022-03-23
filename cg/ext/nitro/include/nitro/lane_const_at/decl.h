#pragma once
#include "../indexed.h"
namespace nitro {

template <typename T, size_t N> struct lane_const_at_expr_impl;

template <typename T, size_t N>
using lane_const_at_expr = typename lane_const_at_expr_impl<T, N>::type;
} // namespace nitro