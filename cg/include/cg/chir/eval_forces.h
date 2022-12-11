#pragma once
#include "chiral_quad.h"
#include <cg/simul/sched.h>
#include <cg/vect/vect.h>

namespace cg::chir {
class eval_forces : public simul::iter_divisible_mixin<eval_forces> {
public:
  real e_chi;

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  vect::const_view<chiral_quad> quads;
  real *V;

public:
  void iter(int idx) const;
  int size() const;
};
} // namespace cg::chir