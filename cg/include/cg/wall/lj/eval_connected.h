#pragma once
#include "data.h"

namespace cg::wall::lj {
class eval_connected {
public:
  real min_dist, breaking_dist, saturation_diff;
  class lj force;

public:
  vect::view<wall> walls;
  vect::set<connection> *conns;
  vect::const_view<vec3r> r;
  vect::view<vec3r> F, wall_F;
  real *V;
  vect::view<bool> is_connected;
  vect::vector<int> *removed;

public:
  void operator()() const;
  void omp_async() const;
  void iter(int conn_idx) const;
};
} // namespace cg::wall::lj