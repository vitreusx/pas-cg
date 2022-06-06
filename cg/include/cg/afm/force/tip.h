#pragma once
#include <cg/types/amp.h>

namespace cg::afm::force {
template <typename E> struct tip_expr { EXPR(res_idx, pull_force) };

class tip : public tip_expr<tip> {
public:
  INST(tip, FIELD(int, res_idx), FIELD(vec3r, pull_force));
};
} // namespace cg::afm::force