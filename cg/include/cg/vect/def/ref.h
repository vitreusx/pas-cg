#pragma once

namespace nitro::def {
template <typename T> using ref = T &;

template <typename T> using const_ref = T const &;
} // namespace nitro::def