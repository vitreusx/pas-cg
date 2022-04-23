#pragma once
#include <cg/vect/vect.h>

namespace cg::qa {
template <typename E> struct cys_neigh_expr : public nitro::ind_expr<E> {
  EXPR_BODY(cys_idx, neigh_idx)
};

template <typename E> struct cys_neigh_auto_expr : public cys_neigh_expr<E> {
  AUTO_EXPR_BODY(cys_idx, neigh_idx)
};

using cys_neigh_base = nitro::tuple_wrapper<int, int>;

class cys_neigh : public cys_neigh_auto_expr<cys_neigh>, public cys_neigh_base {
public:
  using Base = cys_neigh_base;
  using Base::Base;
  using Base::get;
};
} // namespace cg::qa

namespace nitro {
template <>
struct is_indexed_impl<cg::qa::cys_neigh> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::qa::cys_neigh> {
  using type = cg::qa::cys_neigh_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::qa::cys_neigh> {
  using type = cg::qa::cys_neigh_auto_expr<E>;
};
} // namespace nitro