#pragma once
#include "../indexed.h"
namespace nitro {

template <typename T> struct allocator_impl;

template <typename T> using allocator = typename allocator_impl<T>::type;
} // namespace nitro