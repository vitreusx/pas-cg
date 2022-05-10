#pragma once
#include <cg/types/amp.h>

namespace cg::qa {
template <typename E> struct free_pair_expr : public nitro::ind_expr<E> {
  EXPR_BODY(i1, i2, orig_dist)
};

template <typename E> struct free_pair_auto_expr : public free_pair_expr<E> {
  AUTO_EXPR_BODY(i1, i2, orig_dist)
};

using free_pair_base = nitro::tuple_wrapper<int, int, real>;

class free_pair : public free_pair_auto_expr<free_pair>, public free_pair_base {
public:
  using Base = free_pair_base;
  using Base::Base;
  using Base::get;
};
} // namespace cg::qa

namespace nitro {
template <>
struct is_indexed_impl<cg::qa::free_pair> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::qa::free_pair> {
  using type = cg::qa::free_pair_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::qa::free_pair> {
  using type = cg::qa::free_pair_auto_expr<E>;
};
}