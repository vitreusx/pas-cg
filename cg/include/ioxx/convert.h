#pragma once
#include <experimental/type_traits>
#include <sstream>
#include <string>

namespace ioxx {
template <typename V, typename U> struct convert_impl {
  V operator()(U const &x) const {
    auto repr = convert_impl<std::string, U>()(x);
    return convert_impl<V, std::string>()(repr);
  }
};

template <> struct convert_impl<std::string, std::string> {
  std::string operator()(std::string const &s) const;
};

template <typename U> struct convert_impl<std::string, U> {
  std::string operator()(U const &x) const {
    std::stringstream ss;
    ss << x;
    return ss.str();
  }
};

template <typename V> struct convert_impl<V, std::string> {
  V operator()(std::string const &repr) const {
    std::stringstream ss;
    ss << repr;
    V value;
    ss >> value;
    return value;
  }
};

template <typename V, typename U> V convert(U const &x) {
  return convert_impl<V, U>()(x);
}

template <typename V, typename U>
using has_custom_conv_det =
    decltype(std::declval<convert_impl<V, U>>()(std::declval<U const &>()));

template <typename V, typename U>
inline constexpr bool has_custom_conv =
    std::experimental::is_detected_exact_v<V, has_custom_conv_det, U, V>;

template <typename U>
using can_conv_to_str_det =
    decltype(std::declval<std::stringstream>() << std::declval<U const &>());

template <typename U>
inline constexpr bool can_conv_to_str =
    std::experimental::is_detected_v<can_conv_to_str_det, U>;

template <typename V>
using can_conv_from_str_det =
    decltype(std::declval<std::stringstream>() >> std::declval<V &>());

template <typename V>
inline constexpr bool can_conv_from_str =
    std::experimental::is_detected_v<can_conv_from_str_det, V>;

template <typename U, typename V>
inline constexpr bool can_conv = has_custom_conv<U, V> ||
                                 (can_conv_to_str<U> && can_conv_from_str<V>);

} // namespace ioxx
