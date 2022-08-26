#pragma once
#include "data.h"
#include <cg/base_forces/shifted_lj.h>
#include <cg/simul/runtime.h>

namespace cg::wall::lj {
class sift_free : public simul::sliceable_task {
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
  void operator()() const;
  void omp_async() const;
  void iter(int res_idx) const;

  void for_slice(int from, int to) const override;
  int total_size() const override;
  int slice_size() const override;
};
} // namespace cg::wall::lj