#pragma once
#include "data.h"
#include <cg/simul/runtime.h>

namespace cg::wall::harmonic {
class eval_free : public simul::sliceable_task {
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

  void for_slice(int from, int to) const override;
  int total_size() const override;
  int slice_size() const override;
};
} // namespace cg::wall::harmonic