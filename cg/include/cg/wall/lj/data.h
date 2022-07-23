#pragma once
#include <cg/base_forces/lj.h>
#include <cg/types/amp.h>
#include <cg/types/plane.h>

namespace cg::wall::lj {
template <typename E> struct connection_expr {
  EXPR(res_idx, wall_idx, bead_offset, saturation)
};

class connection : public connection_expr<connection> {
public:
  INST(connection, FIELD(int, res_idx), FIELD(int, wall_idx),
       FIELD(vec3r, bead_offset), FIELD(real, saturation))
};

template <typename E> struct candidate_expr { EXPR(res_idx, wall_idx, dist) };

class candidate : public candidate_expr<candidate> {
public:
  INST(candidate, FIELD(int, res_idx), FIELD(int, wall_idx), FIELD(real, dist));
};

struct wall {
  inline wall(plane<real> const &plane, int limit)
      : plane{plane}, num_conns{0}, limit{limit} {}
  
  plane<real> plane;
  int num_conns, limit;
};
} // namespace cg::wall::lj