#pragma once
#include "pair.h"
#include <cg/simul/sched.h>

namespace cg::local_rep {
class eval_forces : public simul::iter_divisible_mixin<eval_forces> {
public:
  real depth, cutoff;

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  vect::const_view<pair> pairs;
  real *V;

public:
  void iter(int idx) const;
  int size() const;
};
} // namespace cg::local_rep