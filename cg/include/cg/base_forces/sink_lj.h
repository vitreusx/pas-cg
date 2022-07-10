#pragma once
#include "lj.h"
#include <cg/types/amp.h>
#include <cg/utils/hash.h>
#include <cg/utils/math.h>
#include <cg/vect/vect.h>

namespace cg {
template <typename E> struct sink_lj_expr {
  EXPR(depth, r_low, r_high)

  auto operator()(real r, real r_inv) const {
    auto r_eff = (r < r_low() ? r_low() : (r < r_high() ? r : r_high()));
    auto s = r_inv * r_eff, s6 = ipow<6>(s), s12 = s6 * s6;
    auto V = depth() * (s12 - 2.0f * s6);
    auto dV_dr = 12.0f * depth() * r_inv * (s6 - s12);
    dV_dr = (r_low() < r && r < r_low()) ? 0 : dV_dr;
    V = clamp<real>(V, -1.0e3, 1.0e3);
    dV_dr = clamp<real>(dV_dr, -1.0e3, 1.0e3);
    return std::make_tuple(V, dV_dr);
  }

  auto cutoff() const {
    return (real)2.0 * r_high();
  }
};

class sink_lj : public sink_lj_expr<sink_lj> {
public:
  INST(sink_lj, FIELD(real, depth), FIELD(real, r_low), FIELD(real, r_high))

  sink_lj() : sink_lj(0.0, 0.0, 0.0){};

  template <typename E>
  explicit sink_lj(lj_expr<E> const &lj)
      : sink_lj(lj.depth(), lj.r_min(), lj.r_min()) {}

  static inline auto cutoff_(real r_high) {
    return (real)2.0 * r_high;
  }
};

struct sink_lj_specs {
  std::optional<quantity> depth, r_low, r_high;
  sink_lj_specs() = default;
  sink_lj_specs(lj_specs const &specs);

  void load(ioxx::xyaml::node const &node);
  operator sink_lj() const;
};

struct ss_sink_lj_specs {
  std::unordered_map<std::pair<amino_acid, amino_acid>, sink_lj_specs> ss_specs;
  void load(ioxx::xyaml::node const &node);
};
} // namespace cg
