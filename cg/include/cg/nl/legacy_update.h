#pragma once
#include "data.h"
#include "exclusion.h"

namespace cg::nl {
class legacy_update {
public:
  real pad_factor;

public:
  nitro::const_view<vec3r> r;
  box<real> const *simul_box;
  data *nl_data;
  int num_particles;
  real const *max_cutoff, *t;
  bool *invalid;
  nitro::const_view<int> chain_idx, seq_idx;
  nitro::const_view<exclusion> all_nat_cont;

public:
  void operator()() const;
};
} // namespace cg::nl