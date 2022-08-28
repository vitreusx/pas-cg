#pragma once
#include "cys_neigh.h"
#include <cg/nl/data.h>
#include <cg/sbox/pbc.h>
#include <cg/simul/runtime.h>
#include <cg/types/amp.h>

namespace cg::qa {
class update_cys_neigh {
public:
  real neigh_radius;

public:
  vect::const_view<vec3r> r;
  sbox::pbc<real> const *simul_box;
  nl::data *nl;
  vect::vector<cys_neigh> *neigh;

public:
  void operator()() const;
};

class reset_cys_neigh {
public:
  vect::const_view<int> cys_indices;
  vect::view<int> neigh_count;

public:
  void operator()() const;
};

class count_cys_neigh : public simul::sliceable_task {
public:
  real neigh_radius;
  vect::const_view<vec3r> r;
  vect::view<int> neigh_count;
  vect::const_view<int> cys_indices;
  vect::vector<cys_neigh> const *neigh;
  sbox::pbc<real> const *simul_box;

public:
  void iter(int idx) const;
  void operator()() const;
  void omp_async() const;

  void for_slice(int from, int to) const override;
  int total_size() const override;

};
} // namespace cg::qa