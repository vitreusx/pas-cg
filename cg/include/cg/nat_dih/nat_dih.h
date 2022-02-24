#pragma once
#include <cg/types/amp.h>
#include <nitro/nitro.h>

namespace cg {
template <typename E> struct nat_dih_expr : public nitro::ind_expr<E> {
  EXPR_BODY(i1, i2, i3, i4, nat_phi)
};

template <typename E> struct nat_dih_auto_expr : public nat_dih_expr<E> {
  AUTO_EXPR_BODY(i1, i2, i3, i4, nat_phi)
};

using nat_dih_base = nitro::tuple_wrapper<int, int, int, int, real>;

class nat_dih : public nat_dih_expr<nat_dih>, public nat_dih_base {
public:
  using Base = nat_dih_base;
  using Base::Base;
  using Base::get;
};
} // namespace cg

namespace nitro {
template <> struct is_indexed_impl<cg::nat_dih> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::nat_dih> {
  using type = cg::nat_dih_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::nat_dih> {
  using type = cg::nat_dih_auto_expr<E>;
};
}; // namespace nitro