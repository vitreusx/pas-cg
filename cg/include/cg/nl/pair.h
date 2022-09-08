#pragma once
#include <cg/types/amp.h>

namespace cg::nl {
template <typename E>
struct pair_expr {
  EXPR(i1, i2, orig_dist, taken)
};

class pair : public pair_expr<pair> {
public:
  INST(pair, FIELD(int, i1), FIELD(int, i2), FIELD(real, orig_dist),
       FIELD(bool, taken))
};
} // namespace cg::nl
