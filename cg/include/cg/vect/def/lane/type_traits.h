#pragma once
#include <type_traits>

namespace nitro::def {
template <typename T>
struct is_lane_like : std::false_type {};

template <typename T>
inline constexpr bool is_lane_like_v = is_lane_like<T>::value;

template <typename T>
struct lane_type;

template <typename T>
using lane_type_t = typename lane_type<T>::type;

template <typename T>
struct lane_size;

template <typename T>
inline constexpr std::size_t lane_size_v = lane_size<T>::value;

template <typename T>
struct lane_width;

template <typename T>
inline constexpr std::size_t
    lane_width_v = 8 * sizeof(lane_type_t<T>) * lane_size_v<T>;

} // namespace nitro::def