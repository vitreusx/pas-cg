#pragma once
#include <experimental/type_traits>

namespace nitro::ind {
template <typename T>
using is_indexed_det = decltype(T::Idxes());

template <typename T>
constexpr bool is_indexed_v =
    std::experimental::is_detected_v<is_indexed_det, T>;

template <typename T>
using subtypes_t = decltype(T::Subtypes());

template <typename T>
using idxes_t = decltype(T::Idxes());

template <typename T, typename E>
using expr_t = decltype(T::template Expr<E>());

template <typename T, typename E>
constexpr bool is_expr_for = std::is_base_of_v<expr_t<T, E>, E>;
} // namespace nitro::ind