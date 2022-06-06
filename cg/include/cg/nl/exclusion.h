#pragma once
#include <cg/types/amp.h>

namespace cg::nl {
template <typename E> struct exclusion_expr { EXPR(i1, i2) };

class exclusion : public exclusion_expr<exclusion> {
public:
  INST(exclusion, FIELD(int, i1), FIELD(int, i2))
};
} // namespace cg::nl
