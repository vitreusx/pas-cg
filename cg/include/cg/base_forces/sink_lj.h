#pragma once
#include "lj.h"
#include <cg/types/amp.h>
#include <cg/utils/math.h>
#include <cg/vect/vect.h>

namespace cg {
template <typename E> struct sink_lj_expr : public nitro::ind_expr<E> {
  EXPR_BODY(depth, r_min, r_max)

  std::tuple<real, real> operator()(real r, real r_inv) const {
    auto r_eff = (r < r_min() ? r_min() : (r < r_max() ? r : r_max()));
    auto s = r_inv * r_eff, s6 = ipow<6>(s), s12 = s6 * s6;
    auto V = depth() * (s12 - 2.0f * s6);
    auto dV_dr = 12.0f * depth() * r_inv * (s6 - s12);
    dV_dr = (r_min() < r && r < r_max()) ? 0 : dV_dr;
    V = clamp<real>(V, -1.0e3, 1.0e3);
    dV_dr = clamp<real>(dV_dr, -1.0e3, 1.0e3);
    return std::make_tuple(V, dV_dr);
  }

  real cutoff() const { return (real)2.0 * r_max(); }
};

template <typename E> struct sink_lj_auto_expr : public sink_lj_expr<E> {
  AUTO_EXPR_BODY(depth, r_min, r_max)
};

using sink_lj_base = nitro::tuple_wrapper<real, real, real>;

class sink_lj : public sink_lj_auto_expr<sink_lj>, public sink_lj_base {
public:
  using Base = sink_lj_base;
  using Base::Base;
  using Base::get;

  sink_lj() : sink_lj(0.0, 0.0, 0.0){};

  template <typename E>
  explicit sink_lj(lj_expr<E> const &lj)
      : sink_lj(lj.depth(), lj.r_min(), lj.r_min()) {}

  static inline real compute_cutoff(real r_max) { return (real)2.0 * r_max; }
};
} // namespace cg

namespace nitro {
template <> struct is_indexed_impl<cg::sink_lj> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::sink_lj> {
  using type = cg::sink_lj_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::sink_lj> {
  using type = cg::sink_lj_auto_expr<E>;
};
}; // namespace nitro