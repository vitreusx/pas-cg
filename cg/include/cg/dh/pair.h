#pragma once
#include <cg/types/amp.h>

namespace cg::dh {
template <typename E> struct pair_expr : public nitro::ind_expr<E> {
  EXPR_BODY(i1, i2, q1_x_q2)
};

template <typename E> struct pair_auto_expr : public pair_expr<E> {
  AUTO_EXPR_BODY(i1, i2, q1_x_q2)
};

using pair_base = nitro::tuple_wrapper<int, int, real>;

class pair : public pair_expr<pair>, public pair_base {
public:
  using Base = pair_base;
  using Base::Base;
  using Base::get;
};
} // namespace cg::dh

namespace nitro {
template <> struct is_indexed_impl<cg::dh::pair> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::dh::pair> {
  using type = cg::dh::pair_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::dh::pair> {
  using type = cg::dh::pair_auto_expr<E>;
};
}; // namespace nitro