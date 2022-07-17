#pragma once
#include <cg/amino/amino_acid.h>
#include <cg/types/amp.h>
#include <cg/utils/random.h>

namespace cg::lang {
class step_base {
public:
  real gamma_factor, temperature;
  solver_real dt;
  int *step_idx;

  vect::view<vec3r> r, v;
  vect::const_view<vec3r> F;
  vect::const_view<amino_acid> atype;
  vect::const_view<real> mass, mass_inv, mass_rsqrt;

  real *t;
  vect::view<vec3sr> y0, y1, y2, y3, y4, y5;
  solver_real *true_t;
  int num_particles;
  rand_gen *gen;
};
} // namespace cg::lang