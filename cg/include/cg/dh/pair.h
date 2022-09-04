#pragma once
#include <cg/types/amp.h>

namespace cg::dh {
template <typename E>
struct pair_expr {
  EXPR(i1, i2, q1_x_q2)
};

class pair : public pair_expr<pair> {
public:
  INST(pair, FIELD(int, i1), FIELD(int, i2), FIELD(real, q1_x_q2))
};
} // namespace cg::dh