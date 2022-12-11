#pragma once
#include "nat_ang.h"
#include <cg/simul/sched.h>

namespace cg::nat_ang {
class eval_forces : public simul::iter_divisible_mixin<eval_forces> {
public:
  real CBA;

public:
  real *V;
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  vect::const_view<nat_ang> angles;

public:
  void iter(int idx) const;
  int size() const;
};
} // namespace cg::nat_ang