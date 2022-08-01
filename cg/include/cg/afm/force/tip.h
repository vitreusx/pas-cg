#pragma once
#include <cg/types/amp.h>
#include <cg/types/avg.h>

namespace cg::afm::force {
template <typename E> struct tip_expr { EXPR(res_idx, pull_force, avg_vel) };

class tip : public tip_expr<tip> {
public:
  INST(tip, FIELD(int, res_idx), FIELD(vec3r, pull_force),
       FIELD(moving_avg<vec3r, real>, avg_vel));
};
} // namespace cg::afm::force