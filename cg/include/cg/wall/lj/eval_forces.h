#pragma once
#include "data.h"
#include <cg/simul/runtime.h>

namespace cg::wall::lj {

class eval_forces : public simul::sliceable_task {
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

  void for_slice(int from, int to) const override;
  int total_size() const override;

};
} // namespace cg::wall::lj