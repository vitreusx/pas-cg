#pragma once
#include "../data.h"
#include <cg/base_forces/lj.h>

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

struct wall : public gen_wall {
  inline wall(cg::plane<real> const &plane, vec3r *F, real avg_t, int limit)
      : gen_wall{plane, F, avg_t}, num_conns{0}, limit{limit} {};
  
  int num_conns, limit;
};
} // namespace cg::wall::lj