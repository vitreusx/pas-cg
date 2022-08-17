#pragma once
#include "eval_forces.h"

namespace cg::pid {
struct active_counts {
  int intra = 0, intra_bb = 0, intra_ss = 0, inter = 0, inter_bb = 0,
      inter_ss = 0;
};

class count_active {
public:
  eval_forces const *eval;
  vect::const_view<int> chain_idx;

public:
  active_counts operator()() const;
};
} // namespace cg::pid