#pragma once
#include <cg/types/amp.h>
#include <cg/utils/math.h>
#include <cg/vect/vect.h>

namespace cg {
template <typename E> struct lj_expr {
  EXPR(depth, r_min)

  decltype(auto) operator()(real r_inv) const {
    auto s = r_inv * r_min(), s6 = ipow<6>(s), s12 = s6 * s6;
    auto V = depth() * (s12 - (real)2.0 * s6);
    auto dV_dr = (real)12.0 * depth() * r_inv * (s6 - s12);
    V = clamp<real>(V, -1.0e3, 1.0e3);
    dV_dr = clamp<real>(dV_dr, -1.0e3, 1.0e3);
    return std::make_tuple(V, dV_dr);
  }

  decltype(auto) cutoff() const {
    return (real)2.0 * r_min();
  }
};

class lj : public lj_expr<lj> {
public:
  INST(lj, FIELD(real, depth), FIELD(real, r_min))

  lj() : lj((real)0.0, (real)0.0){};

  static inline real compute_cutoff(real r_min) {
    return (real)2.0 * r_min;
  }
};
} // namespace cg
