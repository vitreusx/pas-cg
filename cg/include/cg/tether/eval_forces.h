#pragma once
#include "pair.h"
#include <cg/simul/sched.h>

namespace cg::tether {
class eval_forces : public simul::iter_divisible_mixin<eval_forces> {
public:
  real H1, H2;
  real def_length;

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  vect::const_view<pair> tethers;
  real *V;

public:
  void iter(int idx) const;
  int size() const;
};
} // namespace cg::tether