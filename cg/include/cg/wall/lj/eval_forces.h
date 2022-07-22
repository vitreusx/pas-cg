#pragma once
#include "data.h"

namespace cg::wall::lj {

class eval_forces {
public:
  real min_dist, breaking_dist, saturation_diff;

public:
  vect::const_view<wall> walls;
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  real *V;

public:
  void operator()() const;
  void omp_async() const;
  void iter(int idx) const;
};
} // namespace cg::wall::lj