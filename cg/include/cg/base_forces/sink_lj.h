#pragma once
#include "lj.h"
#include <cg/types/amp.h>
#include <cg/utils/math.h>
#include <cg/vect/vect.h>

namespace cg {
template <typename E> struct sink_lj_expr {
  EXPR(depth, r_min, r_max)

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

  real cutoff() const {
    return (real)2.0 * r_max();
  }
};

class sink_lj : public sink_lj_expr<sink_lj> {
public:
  INST(sink_lj, FIELD(real, depth), FIELD(real, r_min), FIELD(real, r_max))

  sink_lj() : sink_lj(0.0, 0.0, 0.0){};

  template <typename E>
  explicit sink_lj(lj_expr<E> const &lj)
      : sink_lj(lj.depth(), lj.r_min(), lj.r_min()) {}

  static inline real compute_cutoff(real r_max) {
    return (real)2.0 * r_max;
  }
};
} // namespace cg
