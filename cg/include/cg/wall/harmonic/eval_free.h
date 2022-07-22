#pragma once
#include "data.h"

namespace cg::wall::harmonic {
class eval_free {
public:
  real depth, min_dist;

public:
  vect::const_view<bool> is_connected;
  vect::const_view<wall> walls;
  vect::const_view<vec3r> r;
  vect::view<vec3r> F, wall_F;
  real *V;

public:
  void operator()() const;
  void omp_async() const;
  void iter(int idx) const;
};
} // namespace cg::wall::harmonic