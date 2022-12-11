#pragma once
#include <cg/sbox/pbc.h>
#include <cg/simul/sched.h>
#include <cg/types/amp.h>

namespace cg::qa {
class prepare_nh : public simul::iter_divisible_mixin<prepare_nh> {
public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> n, h;
  sbox::pbc<real> const *simul_box;
  vect::const_view<int> prev, next;
  int num_particles;

public:
  void iter(int idx) const;
  int size() const;
};
} // namespace cg::qa