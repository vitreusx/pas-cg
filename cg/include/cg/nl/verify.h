#pragma once
#include "data.h"
#include <cg/simul/runtime.h>

namespace cg::nl {
class verify {
public:
  vect::const_view<vec3r> r;
  sbox::pbc<real> const *simul_box;
  bool *invalid, *first_time;
  int num_particles;
  data const *nl_data;
  real *total_disp;

public:
  void omp_async() const;
};
} // namespace cg::nl