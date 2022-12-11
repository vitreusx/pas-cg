#pragma once
#include "pair.h"
#include <cg/sbox/pbc.h>
#include <cg/simul/sched.h>
#include <optional>

namespace cg::const_dh {
class eval_forces : public simul::iter_divisible_mixin<eval_forces> {
public:
  real screen_dist_inv, V_factor, cutoff;
  void set_V_factor(real permittivity);

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  sbox::pbc<real> const *simul_box;
  vect::vector<dh::pair> const *es_pairs = nullptr;
  real *V;

public:
  void iter(int idx) const;
  int size() const;
};
} // namespace cg::const_dh