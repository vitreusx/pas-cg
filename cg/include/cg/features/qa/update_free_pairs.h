#pragma once
#include "free_pair.h"
#include <cg/features/nl/data.h>
#include <cg/types/box.h>

namespace cg::qa {
class update_free_pairs {
public:
  real max_formation_min_dist;

public:
  nitro::vector<vec3r> const *r;
  box<real> const *box;
  nl::data *nl;
  nitro::set<free_pair> *pairs;

public:
  void operator()() const;
};
} // namespace cg::qa