#pragma once

namespace nitro {
template <typename T> struct const_iterator_impl;

template <typename T>
using const_iterator = typename const_iterator_impl<T>::type;
} // namespace nitro