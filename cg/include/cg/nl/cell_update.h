#pragma once
#include "data.h"
#include <utility>
#include <vector>

namespace cg::nl {
class cell_update {
public:
  real pad_factor;

public:
  nitro::const_view<vec3r> r;
  box<real> const *box;
  data *data;
  int num_particles;
  real const *max_cutoff, *t;
  bool *invalid;
  nitro::const_view<int> chain_idx, seq_idx;
  nitro::const_view<pair> all_nat_cont;

  nitro::view<int> res_cell_idx, reordered_idx;
  nitro::vector<int> *num_res_in_cell, *cell_offset;
  std::vector<std::tuple<int, int, real>> *all_pairs;

public:
  void operator()() const;

private:
  int cell_begin(int cell_idx) const;
  int cell_end(int cell_idx) const;
};
} // namespace cg::nl