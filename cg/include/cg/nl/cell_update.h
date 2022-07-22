#pragma once
#include "data.h"
#include "exclusion.h"
#include <utility>
#include <vector>

namespace cg::nl {
class cell_update {
public:
  real cutoff, pad;

public:
  vect::const_view<vec3r> r;
  sbox::pbc<real> const *simul_box;
  data *nl_data;
  int num_particles;
  real *t;
  bool *invalid;
  vect::const_view<int> chain_idx, seq_idx;
  vect::const_view<exclusion> all_nat_cont;

  vect::view<int> res_cell_idx, reordered_idx;
  vect::vector<int> *num_res_in_cell, *cell_offset;
  vect::vector<pair> *all_pairs;

public:
  void operator()() const;

private:
  int cell_begin(int cell_idx) const;
  int cell_end(int cell_idx) const;
};
} // namespace cg::nl