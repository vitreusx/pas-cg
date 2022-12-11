#pragma once
#include "tip.h"
#include <cg/simul/sched.h>

namespace cg::afm::force {
class eval_forces : public simul::iter_divisible_mixin<eval_forces> {
public:
  vect::view<vec3r> F;
  vect::view<tip> afm_tips;
  vect::const_view<vec3r> v;
  real const *t;

public:
  void iter(int idx) const;
  int size() const;
};
} // namespace cg::afm::force