#pragma once
#include <cg/types/amp.h>
#include <cg/vect/vect.h>

namespace cg {
template <typename E> struct nat_dih_expr {
  EXPR(i1, i2, i3, i4, nat_phi)
};

class nat_dih : public nat_dih_expr<nat_dih> {
public:
  INST(nat_dih, FIELD(int, i1), FIELD(int, i2), FIELD(int, i3), FIELD(int, i4),
       FIELD(real, nat_phi))
};
} // namespace cg