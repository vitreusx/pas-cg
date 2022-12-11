#pragma once
#include "pair.h"
#include <cg/nl/data.h>
#include <cg/simul/sched.h>

namespace cg::pauli {

class eval_forces : public simul::iter_divisible_mixin<eval_forces> {
public:
  real depth, r_excl;

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  sbox::pbc<real> const *simul_box;
  vect::vector<pair> const *pairs;
  real *V;

public:
  void iter(int idx) const;
  int size() const;
};
} // namespace cg::pauli