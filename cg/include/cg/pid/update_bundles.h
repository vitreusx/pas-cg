#pragma once
#include "bundle.h"
#include <cg/amino/amino_acid.h>
#include <cg/nl/data.h>
#include <cg/sbox/pbc.h>

namespace cg::pid {
class update_bundles {
private:
  mutable vect::vector<bundle> end_bundles;

public:
  real cutoff;
  bool include4;

public:
  vect::const_view<vec3r> r;
  vect::const_view<int> prev, next, chain_idx, seq_idx;
  vect::const_view<amino_acid> atype;
  sbox::pbc<real> const *simul_box;
  nl::data *nl;
  vect::vector<bundle> *bundles;
  int *fast_iter_end;

public:
  void operator()() const;
};
} // namespace cg::pid