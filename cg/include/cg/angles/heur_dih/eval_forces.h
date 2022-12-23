#pragma once
#include "heur_dih.h"
#include <cg/simul/sched.h>

namespace cg::heur_dih {
class eval_forces : public simul::iter_divisible_mixin<eval_forces> {
public:
  struct coeffs_t {
    solver_real const_[aa_heur_pair::NUM_TYPES];
    solver_real sin[aa_heur_pair::NUM_TYPES];
    solver_real cos[aa_heur_pair::NUM_TYPES];
    solver_real sin2[aa_heur_pair::NUM_TYPES];
    solver_real cos2[aa_heur_pair::NUM_TYPES];
    solver_real sin_cos[aa_heur_pair::NUM_TYPES];
  };
  coeffs_t coeffs;

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  vect::const_view<heur_dih> dihedrals;
  real *V;

public:
  void iter(int idx) const;
  int size() const;
};
} // namespace cg::heur_dih