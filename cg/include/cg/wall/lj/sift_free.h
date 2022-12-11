#pragma once
#include "data.h"
#include <cg/base_forces/shifted_lj.h>
#include <cg/simul/sched.h>

namespace cg::wall::lj {
class sift_free : public simul::iter_divisible_mixin<sift_free> {
public:
  real min_dist;
  shifted_lj force;

public:
  vect::const_view<wall> walls;
  vect::const_view<bool> is_connected;
  vect::vector<candidate> *candidates;
  vect::vector<int> *removed;
  vect::const_view<vec3r> r;
  vect::view<vec3r> F, wall_F;
  real *V;

public:
  void iter(int res_idx) const;
  int size() const;
};
} // namespace cg::wall::lj