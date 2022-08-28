#pragma once
#include <cg/sbox/pbc.h>
#include <cg/simul/runtime.h>
#include <cg/types/amp.h>

namespace cg::qa {
class prepare_nh : public simul::sliceable_task {
public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> n, h;
  sbox::pbc<real> const *simul_box;
  vect::const_view<int> prev, next;
  int num_particles;

public:
  void iter(int idx) const;
  void operator()() const;
  void omp_async() const;

  void for_slice(int from, int to) const override;
  int total_size() const override;

};
} // namespace cg::qa