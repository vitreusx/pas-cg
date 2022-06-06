#pragma once
#include <cg/types/amp.h>

namespace cg::qa {
template <typename E> struct free_pair_expr { EXPR(i1, i2, orig_dist) };

class free_pair : public free_pair_expr<free_pair> {
public:
  INST(free_pair, FIELD(int, i1), FIELD(int, i2), FIELD(real, orig_dist))
};
} // namespace cg::qa
