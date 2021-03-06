#pragma once
#include "free_pair.h"
#include <cg/nl/data.h>
#include <cg/sbox/pbc.h>

namespace cg::qa {
class update_free_pairs {
public:
  real max_formation_min_dist, cutoff;
  bool include4;

public:
  vect::const_view<vec3r> r;
  sbox::pbc<real> const *simul_box;
  nl::data *nl;
  vect::set<free_pair> *pairs;
  vect::const_view<int> chain_idx, seq_idx;

public:
  void operator()() const;
};
} // namespace cg::qa