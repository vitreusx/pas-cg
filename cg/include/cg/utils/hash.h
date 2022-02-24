#pragma once
#include <functional>
#include <utility>

namespace std {
template <typename T1, typename T2> struct hash<std::pair<T1, T2>> {
  size_t operator()(std::pair<T1, T2> const &pair) const {
    auto h1 = std::hash<T1>()(pair.first);
    auto h2 = std::hash<T2>()(pair.second);
    return h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2);
  }
};
} // namespace std