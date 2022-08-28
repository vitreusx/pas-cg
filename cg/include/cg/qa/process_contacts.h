#pragma once
#include "contact.h"
#include "cys_neigh.h"
#include "free_pair.h"
#include <cg/amino/sync_data.h>
#include <cg/base_forces/disulfide.h>
#include <cg/base_forces/sink_lj.h>
#include <cg/sbox/pbc.h>
#include <cg/simul/runtime.h>

namespace cg::qa {
class process_contacts : public simul::sliceable_task {
public:
  vect::vector<sink_lj> ljs;
  real saturation_speed, dt, factor;
  std::optional<real> fixed_cutoff;

  bool disulfide_special_criteria;
  std::optional<disulfide_force> disulfide;
  real ss_def_dist, ss_dist_max_div;
  int max_neigh_count;

  void set_factor(double breaking_factor);

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  sbox::pbc<real> const *simul_box;
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

  void for_slice(int from, int to) const override;
  int total_size() const override;

};
} // namespace cg::qa