#pragma once
#include <cg/types/amp.h>

namespace cg::nl {
template <typename E> struct pair_expr : public nitro::ind_expr<E> {
  EXPR_BODY(i1, i2, orig_dist)
};

template <typename E> struct pair_auto_expr : public pair_expr<E> {
  AUTO_EXPR_BODY(i1, i2, orig_dist)
};

using pair_base = nitro::tuple_wrapper<int, int, real>;

class pair : public pair_auto_expr<pair>, public pair_base {
public:
  using Base = pair_base;
  using Base::Base;
  using Base::get;
};
} // namespace cg::nl

namespace nitro {
template <> struct is_indexed_impl<cg::nl::pair> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::nl::pair> {
  using type = cg::nl::pair_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::nl::pair> {
  using type = cg::nl::pair_auto_expr<E>;
};
}; // namespace nitro