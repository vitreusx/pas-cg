#pragma once
#include "data.h"
#include <cg/simul/runtime.h>

namespace cg::wall::lj {
class eval_connected : public simul::sliceable_task {
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

  void for_slice(int from, int to) const override;
  int total_size() const override;

};
} // namespace cg::wall::lj