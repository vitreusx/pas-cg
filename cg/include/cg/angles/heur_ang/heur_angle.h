#pragma once
#include "../aa_heur_pair.h"
#include <cg/types/amp.h>

namespace cg::heur_ang {
template <typename E> struct heur_ang_expr { EXPR(i1, i2, i3, type) };

class heur_ang : public heur_ang_expr<heur_ang> {
public:
  INST(heur_ang, FIELD(int, i1), FIELD(int, i2), FIELD(int, i3),
       FIELD(aa_heur_pair, type))
};
} // namespace cg::heur_ang
