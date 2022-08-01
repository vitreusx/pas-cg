#pragma once
#include <cg/types/amp.h>
#include <cg/types/avg.h>
#include <cg/types/plane.h>

namespace cg::wall {
struct gen_wall {
  gen_wall(cg::plane<real> const &plane, vec3r *F, real avg_t)
      : plane{plane}, shift{vec3r::Zero()}, F{F}, work{(real)0.0}, avg_F{
                                                                       avg_t} {}

  cg::plane<real> plane;
  vec3r shift, *F;
  real work;
  moving_avg<vec3r, real> avg_F;
};
} // namespace cg::wall