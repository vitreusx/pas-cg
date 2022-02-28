#pragma once
#include "bundle.h"
#include <cg/amino/amino_acid.h>
#include <cg/nl/data.h>
#include <cg/types/box.h>

namespace cg::pid {
class update_bundles {
public:
  real cutoff;

public:
  nitro::const_view<vec3r> r;
  nitro::const_view<int> prev, next;
  nitro::const_view<amino_acid> atype;
  box<real> const *simul_box;
  nl::data const *nl;
  nitro::vector<bundle> *bundles;

public:
  void operator()() const;
};
} // namespace cg::pid