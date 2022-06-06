#pragma once
#include "../aa_heur_pair.h"
#include <cg/types/amp.h>

namespace cg::heur_dih {
template <typename E> struct heur_dih_expr { EXPR(i1, i2, i3, i4, type) };

class heur_dih : public heur_dih_expr<heur_dih> {
public:
  INST(heur_dih, FIELD(int, i1), FIELD(int, i2), FIELD(int, i3), FIELD(int, i4),
       FIELD(aa_heur_pair, type));
};
} // namespace cg::heur_dih
