#pragma once
#include <cg/amino/amino_acid.h>
#include <cg/types/amp.h>
#include <cg/utils/random.h>

namespace cg::lang {
class step {
public:
  static constexpr solver_real f02 = 3.0f / 16.0f, f12 = 251.0f / 360.0f,
                               f32 = 11.0f / 18.0f, f42 = 1.0f / 6.0f,
                               f52 = 1.0f / 60.0f;

  real gamma2, const2;
  solver_real dt, dt_inv, deltsq;
  void set_params(real gamma_factor, real temperature, solver_real dt);

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

public:
  vect::view<vec3r> noise;

public:
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::lang