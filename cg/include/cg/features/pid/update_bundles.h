#pragma once
#include "bundle.h"
#include <cg/features/nl/data.h>
#include <cg/types/amino_acid.h>
#include <cg/types/box.h>

namespace cg::pid {
class update_bundles {
public:
  real cutoff;

public:
  nitro::vector<vec3r> const *r;
  nitro::vector<int> const *prev, *next;
  nitro::vector<amino_acid> atype;
  box<real> const *box;
  nl::data const *nl;
  nitro::vector<bundle> *bundles;

public:
  void operator()() const;
};
} // namespace cg::pid