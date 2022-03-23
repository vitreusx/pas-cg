#pragma once
#include <sstream>
#include <string>

namespace ioxx {
template <typename V, typename U> struct convert_impl {
  V operator()(U const &x) const {
    auto repr = convert_impl<std::string, U>()(x);
    return convert_impl<V, std::string>()(repr);
  }
};

template<> struct convert_impl<std::string, std::string> {
  std::string operator()(std::string const& s) const;
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
} // namespace ioxx
