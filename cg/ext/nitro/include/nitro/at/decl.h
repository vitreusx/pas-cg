#pragma once
#include "../indexed.h"
namespace nitro {

template <typename T> struct at_expr_impl;

template <typename T> using at_expr = typename at_expr_impl<T>::type;
} // namespace nitro