#pragma once
#include <vector>

namespace nitro::debug {
template <typename T>
using iterator = decltype(std::declval<std::vector<T> &>().begin());
}