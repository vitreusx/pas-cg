#pragma once
#include <cg/types/amp.h>
#include <cg/types/avg.h>
#include <cg/types/plane.h>

namespace cg::wall {
class log_forces {
public:
  vec3r const *neg_z_force, *pos_z_force;
  plane<real> const *neg_z_plane, *pos_z_plane;
  real const *t;
  moving_avg<real, real> *avg_z_force;

public:
  void operator()() const;
};
} // namespace cg::wall