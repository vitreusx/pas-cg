#pragma once
#include "candidate.h"
#include "free_pair.h"
#include <cg/amino/amino_acid.h>
#include <cg/amino/sync_data.h>
#include <cg/types/box.h>

namespace cg::qa {
class loop_over_candidates {
public:
  real min_abs_cos_hr, min_abs_cos_hh, max_cos_nr;
  real req_min_dist[contact_type::NUM_TYPES];
  real max_req_dist, rep_cutoff, rep_depth;
  polarization_type ptype[amino_acid::NUM_TYPES];
  real formation_tolerance;

  bool disulfide_special_criteria;
  int max_neigh_count;

public:
  vect::const_view<vec3r> r, n, h;
  vect::view<vec3r> F;
  box<real> const *simul_box;
  vect::const_view<amino_acid> atype;
  vect::const_view<sync_data> sync;
  real *V, *total_disp;

  vect::set<free_pair> const *free_pairs;
  vect::vector<candidate> *candidates;

  vect::const_view<int> neigh;
  vect::const_view<bool> part_of_ssbond;

public:
  void iter(int idx) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::qa