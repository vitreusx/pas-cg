#pragma once
#include <cg/types/amp.h>
#include <nitro/nitro.h>

namespace cg::nat_cont {
template <typename E> struct nat_cont_expr : public nitro::ind_expr<E> {
  EXPR_BODY(i1, i2, nat_dist)
};

template <typename E> struct nat_cont_auto_expr : public nat_cont_expr<E> {
  AUTO_EXPR_BODY(i1, i2, nat_dist)
};

using nat_cont_base = nitro::tuple_wrapper<int, int, real>;

class nat_cont : public nat_cont_expr<nat_cont>, public nat_cont_base {
public:
  using Base = nat_cont_base;
  using Base::Base;
  using Base::get;
};
} // namespace cg::nat_cont

namespace nitro {
template <>
struct is_indexed_impl<cg::nat_cont::nat_cont> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::nat_cont::nat_cont> {
  using type = cg::nat_cont::nat_cont_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::nat_cont::nat_cont> {
  using type = cg::nat_cont::nat_cont_auto_expr<E>;
};
}; // namespace nitro