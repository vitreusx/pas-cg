#pragma once
#include "contact.h"
#include "cys_neigh.h"
#include "free_pair.h"
#include <cg/amino/sync_data.h>
#include <cg/base_forces/disulfide.h>
#include <cg/base_forces/sink_lj.h>
#include <cg/types/box.h>

namespace cg::qa {
class process_contacts {
public:
  vect::vector<sink_lj> ljs;
  real cycle_time, cycle_time_inv;
  real factor;
  std::optional<real> fixed_cutoff;

  bool disulfide_special_criteria;
  std::optional<disulfide_force> disulfide;
  real ss_def_dist, ss_dist_max_div;
  int max_neigh_count;

  void set_factor(double breaking_factor);
  real saturation_value(contact const &cont) const;

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  box<real> const *simul_box;
  vect::set<contact> *contacts;
  real *V, *t;
  vect::view<sync_data> sync;
  vect::set<free_pair> *free_pairs;
  vect::const_view<int> neigh;
  vect::view<bool> part_of_ssbond;
  vect::vector<int> *removed;

public:
  void iter(int idx) const;

  void operator()() const;
  void omp_async() const;
};
} // namespace cg::qa