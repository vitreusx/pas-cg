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
  nitro::const_view<vec3r> r;
  box<real> const *simul_box;
  nl::data *nl;
  nitro::vector<cys_neigh> *neigh;

public:
  void operator()() const;
};

class count_cys_neigh {
public:
  real neigh_radius;
  nitro::const_view<vec3r> r;
  nitro::view<int> neigh_count;
  nitro::const_view<int> cys_indices;
  nitro::vector<cys_neigh> const *neigh;
  box<real> const *simul_box;

public:
  void operator()() const;
};
} // namespace cg::qa