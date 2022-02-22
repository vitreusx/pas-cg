#pragma once
#include "pair.h"
#include <cg/features/nl/data.h>
#include <cg/types/amino_acid.h>
#include <cg/types/amp.h>
#include <cg/types/box.h>

namespace cg::dh {
class update_pairs {
public:
  real cutoff;
  real q[amino_acid::NUM_TYPES];

public:
  nitro::vector<vec3r> const *r;
  box<real> const *box;
  nl::data const *nl;
  nitro::vector<pair> *pairs;
  nitro::vector<amino_acid> const *atype;

public:
  void operator()() const;
};
} // namespace cg::dh