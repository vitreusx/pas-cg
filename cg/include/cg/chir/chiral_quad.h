#pragma once
#include <cg/types/amp.h>

namespace cg::chir {
template <typename E> struct chiral_quad_expr : public nitro::ind_expr<E> {
  EXPR_BODY(i1, i2, i3, i4, nat_chir, nat_factor)
};

template <typename E>
struct chiral_quad_auto_expr : public chiral_quad_expr<E> {
  AUTO_EXPR_BODY(i1, i2, i3, i4, nat_chir, nat_factor)
};

using chiral_quad_base = nitro::tuple_wrapper<int, int, int, int, real, real>;

class chiral_quad : public chiral_quad_expr<chiral_quad>,
                    public chiral_quad_base {
public:
  using Base = chiral_quad_base;
  using Base::Base;
  using Base::get;
};
} // namespace cg::chir

namespace nitro {
template <>
struct is_indexed_impl<cg::chir::chiral_quad> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::chir::chiral_quad> {
  using type = cg::chir::chiral_quad_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::chir::chiral_quad> {
  using type = cg::chir::chiral_quad_auto_expr<E>;
};
}