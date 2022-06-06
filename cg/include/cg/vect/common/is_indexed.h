#pragma once
#include <experimental/type_traits>

template <typename T> using is_indexed_det = decltype(T::Idxes());

template <typename T>
constexpr bool is_indexed_v =
    std::experimental::is_detected_v<is_indexed_det, T>;