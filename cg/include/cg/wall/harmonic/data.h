#pragma once
#include <cg/types/amp.h>
#include <cg/types/plane.h>

namespace cg::wall::harmonic {
template <typename E> struct connection_expr {
  EXPR(wall_idx, res_idx, bead_offset, nat_dist)
};

class connection : public connection_expr<connection> {
public:
  INST(connection, FIELD(int, wall_idx), FIELD(int, res_idx),
       FIELD(vec3r, bead_offset), FIELD(real, nat_dist))
};

struct wall {
  inline wall(cg::plane<real> const &plane, int limit)
      : plane{plane}, limit{limit} {}

  cg::plane<real> plane;
  int limit;
};
} // namespace cg::wall::harmonic