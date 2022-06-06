#pragma once
#include <deque>

namespace nitro::debug {
template <typename T>
using const_iterator = decltype(std::declval<std::deque<T> const &>().begin());
}