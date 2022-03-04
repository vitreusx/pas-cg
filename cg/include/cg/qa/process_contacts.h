#pragma once
#include "contact.h"
#include "cys_neigh.h"
#include "free_pair.h"
#include "lj_variants.h"
#include <cg/amino/sync_data.h>
#include <cg/base_forces/disulfide.h>
#include <cg/types/box.h>

namespace cg::qa {
class process_contacts {
public:
  lj_variants ljs;
  real cycle_time, cycle_time_inv;
  real factor;

  bool disulfide_special_criteria;
  std::optional<disulfide_force> disulfide;
  real ss_def_dist, ss_dist_max_div;
  int max_neigh_count;

  void set_factor(double breaking_factor);
  real saturation_value(contact const &cont) const;

public:
  nitro::const_view<vec3r> r;
  nitro::view<vec3r> F;
  box<real> const *simul_box;
  nitro::set<contact> *contacts;
  real *V, *t;
  nitro::view<sync_data> sync;
  nitro::set<free_pair> *free_pairs;
  nitro::const_view<int> neigh;

public:
  void iter(int idx) const;

  void operator()() const;
  void omp_async() const;
};
} // namespace cg::qa