#pragma once
#include "data.h"
#include "exclusion.h"
#include <cg/simul/runtime.h>

namespace cg::nl {
class legacy_update {
public:
  real cutoff, pad;

public:
  vect::const_view<vec3r> r;
  sbox::pbc<real> const *simul_box;
  data *nl_data;
  real const *t;
  bool *invalid;
  vect::const_view<int> chain_idx, seq_idx, idxes;

public:
  void omp_async() const;
};
} // namespace cg::nl