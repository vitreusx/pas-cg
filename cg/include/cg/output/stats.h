#pragma once
#include "hook.h"
#include <cg/amino/amino_acid.h>
#include <cg/types/amp.h>
#include <cg/utils/quantity.h>

namespace cg::out {
struct gyration_stuff {
  vec3r lambda;
  quantity radius_of_gyration() const;
  quantity asphericity() const;
};

class add_stats : public hook {
public:
  nitro::const_view<vec3r> orig_r, r, v;
  nitro::const_view<real> mass;
  nitro::const_view<amino_acid> atype;
  nitro::const_view<int> chain_first, chain_last, chain_idx;
  real *V, *t;

public:
  gyration_stuff gyration(nitro::vector<int> const &indices) const;
  real rmsd(nitro::vector<int> const &indices) const;

  void report_to(report_data &report) const override;
};
} // namespace cg::out