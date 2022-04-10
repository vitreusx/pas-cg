#pragma once

namespace nitro {
template <typename T> struct iterator_impl;

template <typename T> using iterator = typename iterator_impl<T>::type;
} // namespace nitro