#pragma once
#include <cg/amino/amino_acid.h>
#include <cg/types/amp.h>
#include <cg/utils/random.h>

namespace cg::lang {

class step {
public:
  real gamma_factor, temperature;
  solver_real dt;

public:
  nitro::view<vec3r> r, v;
  nitro::const_view<vec3r> F;
  nitro::const_view<amino_acid> atype;
  nitro::const_view<real> mass, mass_inv, mass_rsqrt;

  real *t;
  nitro::view<vec3sr> y0, y1, y2, y3, y4, y5;
  solver_real *true_t;
  int num_particles;
  rand_gen *gen;

public:
  void operator()() const;
};
} // namespace cg::lang