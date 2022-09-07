#pragma once
#include "type_traits.h"
#include <array>

namespace nitro::def {
template <typename T, std::size_t N>
using array_lane = std::array<T, N>;

template <typename T, std::size_t N>
struct is_lane_like<std::array<T, N>> : std::true_type {};

template <typename T, std::size_t N>
struct lane_type<std::array<T, N>> {
  using type = T;
};

template <typename T, std::size_t N>
struct lane_size<std::array<T, N>> {
  static constexpr std::size_t value = N;
};

} // namespace nitro::def