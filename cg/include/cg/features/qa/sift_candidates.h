#pragma once
#include "candidate.h"
#include "free_pair.h"
#include "sync_data.h"
#include <cg/types/amino_acid.h>
#include <cg/types/box.h>

namespace cg::qa {
class sift_candidates {
public:
  real min_abs_cos_hr, min_abs_cos_hh, max_cos_nr;
  real req_min_dist[contact_type::NUM_TYPES];
  polarization_type ptype[amino_acid::NUM_TYPES];

public:
  nitro::vector<vec3r> const *r, *n, *h;
  box<real> const *box;
  nitro::vector<amino_acid> const *atype;
  nitro::vector<sync_data> const *sync;

  nitro::set<free_pair> const *free_pairs;
  nitro::vector<candidate> *candidates;

public:
  void iter(int idx) const;
  void operator()() const;

  void omp_prep() const;
  void omp_async() const;
};
} // namespace cg::qa