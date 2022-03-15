#pragma once
#include "hook.h"
#include <cg/amino/amino_acid.h>
#include <cg/types/amp.h>

namespace cg::out {
class add_stats : public hook {
public:
  nitro::const_view<vec3r> r, v;
  nitro::const_view<real> mass;
  nitro::const_view<amino_acid> atype;
  nitro::const_view<int> chain_first, chain_last;
  real *V, *t;

public:
  void report_to(report_data &report) const override;
};
} // namespace cg::out