#pragma once
#include <cg/types/amp.h>

namespace cg::local_rep {
template <typename E> struct pair_expr { EXPR(i1, i2) };

class pair : public pair_expr<pair> {
public:
  INST(pair, FIELD(int, i1), FIELD(int, i2))
};
} // namespace cg::local_rep