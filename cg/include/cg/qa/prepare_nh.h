#pragma once
#include <cg/types/amp.h>
#include <cg/types/box.h>

namespace cg::qa {
class prepare_nh {
public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> n, h;
  box<real> const *simul_box;
  vect::const_view<int> prev, next;
  int num_particles;

public:
  void iter(int idx) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::qa