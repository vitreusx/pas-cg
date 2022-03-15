#pragma once
#include "data.h"

namespace cg::nl {
class verify {
public:
  nitro::const_view<vec3r> r;
  box<real> const *simul_box;
  bool *invalid, *first_time;
  int num_particles;
  data const *nl_data;
  real *total_disp;

public:
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::nl