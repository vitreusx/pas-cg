#pragma once
#include "pair.h"
#include <cg/nl/data.h>
#include <cg/types/box.h>

namespace cg::pauli {
class update_pairs {
public:
  real r_excl;

public:
  vect::const_view<vec3r> r;
  box<real> const *simul_box;
  nl::data const *nl;
  vect::vector<pair> *pairs;

public:
  void operator()() const;
};
} // namespace cg::pauli