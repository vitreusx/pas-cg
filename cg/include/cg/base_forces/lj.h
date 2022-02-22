#pragma once
#include <cg/types/amp.h>
#include <cg/utils/math.h>
#include <nitro/nitro.h>

namespace cg {
template <typename E> struct lj_expr : public nitro::ind_expr<E> {
  EXPR_BODY(depth, r_min)

  decltype(auto) operator()(real r_inv) const {
    auto s = r_inv * r_min(), s6 = ipow<6>(s), s12 = s6 * s6;
    auto V = depth() * (s12 - (real)2.0 * s6);
    auto dV_dr = (real)12.0 * depth() * r_inv * (s6 - s12);
    V = clamp<real>(V, -1.0e3, 1.0e3);
    dV_dr = clamp<real>(dV_dr, -1.0e3, 1.0e3);
    return std::make_tuple(V, dV_dr);
  }

  decltype(auto) cutoff() const { return (real)2.0 * r_min(); }
};

template <typename E> struct lj_auto_expr : public lj_expr<E> {
  AUTO_EXPR_BODY(depth, r_min)
};

using lj_base = nitro::tuple_wrapper<int, real>;

class lj : public lj_auto_expr<lj>, public lj_base {
public:
  using Base = lj_base;
  using Base::Base;
  using Base::get;

  static inline real cutoff(real r_min) { return (real)2.0 * r_min; }
};
} // namespace cg

namespace nitro {
template <> struct is_indexed_impl<cg::lj> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::lj> {
  using type = cg::lj_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::lj> {
  using type = cg::lj_auto_expr<E>;
};
}; // namespace nitro