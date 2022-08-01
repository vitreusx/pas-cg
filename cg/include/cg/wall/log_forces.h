#pragma once
#include "data.h"
#include <cg/types/amp.h>
#include <cg/types/avg.h>
#include <cg/types/plane.h>

namespace cg::wall {
class log_forces {
public:
  real const *t;
  vect::view<wall::gen_wall *> walls;
  wall::gen_wall *neg_z, *pos_z;
  moving_avg<real, real> *avg_z_force;

public:
  void operator()() const;
};
} // namespace cg::wall