#pragma once
#include "data.h"
#include <cg/simul/sched.h>

namespace cg::wall::harmonic {
class eval_connected : public simul::iter_divisible_mixin<eval_connected> {
public:
  real HH1;

public:
  vect::const_view<wall> walls;
  vect::const_view<connection> conns;
  vect::const_view<vec3r> r;
  vect::view<vec3r> F, wall_F;
  real *V;

public:
  void iter(int idx) const;
  int size() const;
};
} // namespace cg::wall::harmonic