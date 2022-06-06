#pragma once
#include "cys_neigh.h"
#include <cg/nl/data.h>
#include <cg/types/amp.h>
#include <cg/types/box.h>

namespace cg::qa {
class update_cys_neigh {
public:
  real neigh_radius;

public:
  vect::const_view<vec3r> r;
  box<real> const *simul_box;
  nl::data *nl;
  vect::vector<cys_neigh> *neigh;

public:
  void operator()() const;
};

class count_cys_neigh {
public:
  real neigh_radius;
  vect::const_view<vec3r> r;
  vect::view<int> neigh_count;
  vect::const_view<int> cys_indices;
  vect::vector<cys_neigh> const *neigh;
  box<real> const *simul_box;

public:
  void operator()() const;
  void omp_reset() const;
  void omp_async() const;
};
} // namespace cg::qa