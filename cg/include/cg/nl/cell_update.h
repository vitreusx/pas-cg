#pragma once
#include "data.h"
#include "exclusion.h"
#include <cg/simul/runtime.h>
#include <unordered_map>

namespace cg::nl {
class cell_update {
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
  vec3r *shared_r_min, *shared_r_max;
  int *total_pairs, *cur_pairs_offset;
  vect::view<std::pair<int, int>> cell_idxes;
  vect::vector<std::tuple<int, int, int>> *unique_cells;
  std::unordered_map<int, int> *unique_cell_idxes;
  mutable vect::vector<std::tuple<int, int, int>> offsets;
  mutable vect::vector<nl::pair> local_pairs;

public:
  void omp_async() const;
};
} // namespace cg::nl