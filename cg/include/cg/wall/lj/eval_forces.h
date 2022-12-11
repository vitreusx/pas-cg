#pragma once
#include "data.h"
#include <cg/simul/sched.h>

namespace cg::wall::lj {

class eval_forces : public simul::iter_divisible_mixin<eval_forces> {
public:
  real min_dist, breaking_dist, saturation_diff;

public:
  vect::const_view<wall> walls;
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  real *V;

public:
  void iter(int idx) const;
  int size() const;
};
} // namespace cg::wall::lj