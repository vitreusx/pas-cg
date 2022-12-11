#pragma once
#include "data.h"
#include <cg/simul/sched.h>

namespace cg::wall::solid {
class eval_forces : public simul::iter_divisible_mixin<eval_forces> {
public:
  real min_dist;

public:
  vect::const_view<wall> walls;
  real depth;

  vect::const_view<vec3r> r;
  vect::view<vec3r> F, wall_F;

public:
  void iter(int idx) const;
  int size() const;
};
} // namespace cg::wall::solid