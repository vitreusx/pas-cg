#pragma once
#include <cg/types/amp.h>

namespace cg::afm::vel {
template <typename E> struct tip_expr {
  EXPR(res_idx, afm_orig, ref_t, afm_vel);
};

class tip : public tip_expr<tip> {
public:
  INST(tip, FIELD(int, res_idx), FIELD(vec3r, afm_orig), FIELD(real, ref_t),
       FIELD(vec3r, afm_vel));
};
} // namespace cg::afm::vel
