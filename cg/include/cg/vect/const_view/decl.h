#pragma once

namespace nitro {
template <typename T, typename Idx> struct const_view_impl;

template <typename T, typename Idx = int>
using const_view = typename const_view_impl<T, Idx>::type;
} // namespace nitro