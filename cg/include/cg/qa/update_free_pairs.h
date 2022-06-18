#pragma once
#include "free_pair.h"
#include <cg/nl/data.h>
#include <cg/types/box.h>

namespace cg::qa {
class update_free_pairs {
public:
  real max_formation_min_dist;
  std::optional<real> fixed_cutoff;
  bool include4;

public:
  vect::const_view<vec3r> r;
  box<real> const *simul_box;
  nl::data *nl;
  vect::set<free_pair> *pairs;
  vect::const_view<int> chain_idx, seq_idx;

public:
  void operator()() const;
};
} // namespace cg::qa