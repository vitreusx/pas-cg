#pragma once
#include <cg/vect/vect.h>

namespace cg::qa {
template <typename E> struct cys_neigh_expr { EXPR(cys_idx, neigh_idx) };

class cys_neigh : public cys_neigh_expr<cys_neigh> {
public:
  INST(cys_neigh, FIELD(int, cys_idx), FIELD(int, neigh_idx))
};
} // namespace cg::qa
