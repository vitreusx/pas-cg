#pragma once
#include "tip.h"
#include <cg/base_forces/harmonic.h>

namespace cg::afm::vel {
class eval_forces {
public:
  harmonic afm_force;

public:
  nitro::const_view<vec3r> r;
  nitro::view<vec3r> F;
  real *t;
  nitro::const_view<tip> afm_tips;

public:
  template <typename E> void iter(tip_expr<E> const &tip) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::afm::vel