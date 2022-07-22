#pragma once
#include <cg/types/amp.h>
#include <cg/types/plane.h>

namespace cg::wall::solid {
class eval_forces {
public:
  real min_dist;

public:
  vect::const_view<plane<real>> walls;
  real depth;

  vect::const_view<vec3r> r;
  vect::view<vec3r> F, wall_F;

public:
  void operator()() const;
  void omp_async() const;
  void iter(int idx) const;
};
} // namespace cg::wall::solid