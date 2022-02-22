#pragma once
#include <cg/types/amp.h>
#include <cg/utils/random.h>

namespace cg::lang {

class step {
public:
  real gamma_factor, temperature;
  solver_real dt;

public:
  nitro::vector<vec3r> *r, *v;
  nitro::vector<vec3r> const *F;
  nitro::vector<real> const *mass, *mass_inv, *mass_rsqrt;

  real *t;
  nitro::vector<vec3sr> *y0, *y1, *y2, *y3, *y4, *y5;
  solver_real *true_t;
  int num_particles;
  rand_gen *gen;

public:
  void operator()() const;
};
} // namespace cg::lang