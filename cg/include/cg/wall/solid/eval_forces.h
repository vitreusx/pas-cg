#pragma once
#include "data.h"
#include <cg/simul/runtime.h>

namespace cg::wall::solid {
class eval_forces : public simul::sliceable_task {
public:
  real min_dist;

public:
  vect::const_view<wall> walls;
  real depth;

  vect::const_view<vec3r> r;
  vect::view<vec3r> F, wall_F;

public:
  void operator()() const;
  void omp_async() const;
  void iter(int idx) const;

  void for_slice(int from, int to) const override;
  int total_size() const override;

};
} // namespace cg::wall::solid