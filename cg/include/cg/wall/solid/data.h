#pragma once
#include "../data.h"

namespace cg::wall::solid {
struct wall : public gen_wall {
  wall(cg::plane<real> const &plane, vec3r *F, real avg_t)
      : gen_wall{plane, F, avg_t} {}
};
} // namespace cg::wall::solid