#pragma once

namespace nitro {
template <typename T, typename Idx> struct view_impl;

template <typename T, typename Idx = int>
using view = typename view_impl<T, Idx>::type;
} // namespace nitro