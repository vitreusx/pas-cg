#pragma once
#include "../data.h"

namespace cg::wall::harmonic {
template <typename E> struct connection_expr {
  EXPR(wall_idx, res_idx, bead_offset, nat_dist)
};

class connection : public connection_expr<connection> {
public:
  INST(connection, FIELD(int, wall_idx), FIELD(int, res_idx),
       FIELD(vec3r, bead_offset), FIELD(real, nat_dist))
};

struct wall : public gen_wall {
  inline wall(cg::plane<real> const &plane, vec3r *F, real avg_t, int limit)
      : gen_wall{plane, F, avg_t}, limit{limit} {}

  int limit;
};
} // namespace cg::wall::harmonic