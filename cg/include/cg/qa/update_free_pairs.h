#pragma once
#include "free_pair.h"
#include <cg/nl/data.h>
#include <cg/types/box.h>

namespace cg::qa {
class update_free_pairs {
public:
  real max_formation_min_dist;
  bool include4;

public:
  nitro::const_view<vec3r> r;
  box<real> const *simul_box;
  nl::data *nl;
  nitro::set<free_pair> *pairs;
  nitro::const_view<int> chain_idx, seq_idx;

public:
  void operator()() const;
};
} // namespace cg::qa