#pragma once
#include "tip.h"
#include <cg/base_forces/harmonic.h>
#include <cg/simul/sched.h>

namespace cg::afm::vel {
class eval_forces : public simul::iter_divisible_mixin<eval_forces> {
public:
  harmonic afm_force;

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  real *t, *V;
  vect::view<tip> afm_tips;

public:
  void iter(int idx) const;
  int size() const;

  real compute_force(vel::tip const &tip) const;
};
} // namespace cg::afm::vel