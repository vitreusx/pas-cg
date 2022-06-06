#pragma once
#include "bundle.h"
#include <cg/amino/amino_acid.h>
#include <cg/nl/data.h>
#include <cg/types/box.h>

namespace cg::pid {
class update_bundles {
public:
  real cutoff;
  bool include4;

public:
  vect::const_view<vec3r> r;
  vect::const_view<int> prev, next, chain_idx, seq_idx;
  vect::const_view<amino_acid> atype;
  box<real> const *simul_box;
  nl::data const *nl;
  vect::vector<bundle> *bundles;

public:
  void operator()() const;
};
} // namespace cg::pid