#pragma once
#include <cg/types/amp.h>

namespace cg::pauli {
template <typename E> struct pair_expr : public nitro::ind_expr<E> {
  EXPR_BODY(i1, i2)
};

template <typename E> struct pair_auto_expr : public pair_expr<E> {
  AUTO_EXPR_BODY(i1, i2)
};

using pair_base = nitro::tuple_wrapper<int, int>;

class pair : public pair_expr<pair>, public pair_base {
public:
  using Base = pair_base;
  using Base::Base;
  using Base::get;
};
} // namespace cg::pauli

namespace nitro {
template <> struct is_indexed_impl<cg::pauli::pair> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::pauli::pair> {
  using type = cg::pauli::pair_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::pauli::pair> {
  using type = cg::pauli::pair_auto_expr<E>;
};
}