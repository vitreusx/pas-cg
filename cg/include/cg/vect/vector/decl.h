#pragma once
#include "../allocator.h"

namespace nitro {
template <typename T, typename Alloc, typename Idx> struct vector_impl;

template <typename T, typename Alloc = allocator<T>, typename Idx = int>
using vector = typename vector_impl<T, Alloc, Idx>::type;
} // namespace nitro
