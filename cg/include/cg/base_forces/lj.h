#pragma once
#include <cg/amino/amino_acid.h>
#include <cg/files/files.h>
#include <cg/types/amp.h>
#include <cg/utils/hash.h>
#include <cg/utils/math.h>
#include <cg/utils/quantity.h>
#include <cg/vect/vect.h>

#define C216 1.122462048309373
#define C216_INV 0.8908987181403393

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

struct lj_specs {
  std::optional<quantity> depth, r_min;
  operator lj() const;
  void load(ioxx::xyaml::node const &node);
};

struct ss_lj_specs {
  std::unordered_map<std::pair<amino_acid, amino_acid>, lj_specs> ss_specs;
  void load(ioxx::xyaml::node const &node);
};
} // namespace cg
