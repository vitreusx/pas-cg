#pragma once
#include "pair.h"
#include <cg/amino/amino_acid.h>
#include <cg/nl/data.h>
#include <cg/types/amp.h>
#include <cg/types/box.h>

namespace cg::dh {
class update_pairs {
public:
  real cutoff;
  real q[amino_acid::NUM_TYPES];

public:
  vect::const_view<vec3r> r;
  box<real> const *simul_box;
  nl::data const *nl;
  vect::vector<pair> *pairs;
  vect::const_view<amino_acid> atype;

public:
  void operator()() const;
};
} // namespace cg::dh