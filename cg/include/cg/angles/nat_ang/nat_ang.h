#pragma once
#include <cg/types/amp.h>
#include <cg/vect/vect.h>

namespace cg::nat_ang {
template <typename E> struct nat_ang_expr { EXPR(i1, i2, i3, nat_theta) };

class nat_ang : public nat_ang_expr<nat_ang> {
public:
  INST(nat_ang, FIELD(int, i1), FIELD(int, i2), FIELD(int, i3),
       FIELD(real, nat_theta))
};
} // namespace cg::nat_ang
