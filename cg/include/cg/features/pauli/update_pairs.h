#pragma once
#include "pair.h"
#include <cg/features/nl/data.h>
#include <cg/types/box.h>

namespace cg::pauli {
class update_pairs {
public:
  real r_excl;

public:
  nitro::vector<vec3r> const *r;
  box<real> const *box;
  nl::data const *nl;
  nitro::vector<pair> *pairs;

public:
  void operator()() const;
};
} // namespace cg::pauli