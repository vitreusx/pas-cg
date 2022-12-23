#pragma once
#include "heur_angle.h"
#include <cg/simul/sched.h>

namespace cg::heur_ang {
class eval_forces : public simul::iter_divisible_mixin<eval_forces> {
public:
  static constexpr int POLY_DEG = 6, NUM_TYPES = 9;
  solver_real poly_coeffs[POLY_DEG + 1][NUM_TYPES];

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  vect::const_view<heur_ang> angles;
  real *V;

public:
  void iter(int idx) const;
  int size() const;
};
} // namespace cg::heur_ang