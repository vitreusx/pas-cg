#pragma once
#include <cg/types/amp.h>

namespace cg::chir {
template <typename E> struct chiral_quad_expr {
  EXPR(i1, i2, i3, i4, nat_chir, nat_factor);
};

class chiral_quad : public chiral_quad_expr<chiral_quad> {
public:
  INST(chiral_quad, FIELD(int, i1), FIELD(int, i2), FIELD(int, i3),
       FIELD(int, i4), FIELD(real, nat_chir), FIELD(real, nat_factor));
};
} // namespace cg::chir
