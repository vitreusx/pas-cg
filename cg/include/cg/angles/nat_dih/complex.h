#pragma once
#include "nat_dih.h"
#include <cg/simul/sched.h>

namespace cg::cnd {
class eval_forces : public simul::iter_divisible_mixin<eval_forces> {
public:
  real CDA, CDB;
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  vect::const_view<nat_dih> dihedrals;
  real *V;

public:
  void iter(int idx) const;
  int size() const;
};
} // namespace cg::cnd
