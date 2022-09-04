#pragma once
#include <cg/types/amp.h>
#include <cg/types/avg.h>

namespace cg::afm::vel {
template <typename E>
struct tip_expr {
  EXPR(res_idx, afm_orig, ref_t, afm_vel, avg_force, avg_perp_force);
};

class tip : public tip_expr<tip> {
public:
  INST(tip, FIELD(int, res_idx), FIELD(vec3r, afm_orig), FIELD(real, ref_t),
       FIELD(vec3r, afm_vel), FIELD(moving_avg<vec3r, real>, avg_force),
       FIELD(moving_avg<real, real>, avg_perp_force))
};
} // namespace cg::afm::vel
