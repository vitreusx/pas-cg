#pragma once
#include <cg/types/amp.h>
#include <cg/vect/vect.h>

namespace cg::nat_ang {
template <typename E> struct nat_ang_expr : public nitro::ind_expr<E> {
  EXPR_BODY(i1, i2, i3, nat_theta)
};

template <typename E> struct nat_ang_auto_expr : public nat_ang_expr<E> {
  AUTO_EXPR_BODY(i1, i2, i3, nat_theta)
};

using nat_ang_base = nitro::tuple_wrapper<int, int, int, real>;

class nat_ang : public nat_ang_expr<nat_ang>, public nat_ang_base {
public:
  using Base = nat_ang_base;
  using Base::Base;
  using Base::get;
};
} // namespace cg::nat_ang

namespace nitro {
template <>
struct is_indexed_impl<cg::nat_ang::nat_ang> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::nat_ang::nat_ang> {
  using type = cg::nat_ang::nat_ang_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::nat_ang::nat_ang> {
  using type = cg::nat_ang::nat_ang_auto_expr<E>;
};
}