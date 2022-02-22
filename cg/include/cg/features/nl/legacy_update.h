#pragma once
#include "data.h"

namespace cg::nl {
class legacy_update {
public:
  real pad_factor;

public:
  nitro::vector<vec3r> const *r;
  box<real> const *box;
  data *data;
  int num_particles;
  real const *max_cutoff, *t;
  bool *invalid;
  nitro::vector<int> const *chain_idx, *seq_idx;
  nitro::vector<pair> const *all_nat_cont;

public:
  void operator()() const;
};
} // namespace cg::nl