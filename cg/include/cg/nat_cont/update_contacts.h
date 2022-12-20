#pragma once
#include "nat_cont.h"
#include <cg/nl/data.h>
#include <cg/sbox/pbc.h>

namespace cg::nat_cont {
class update_contacts {
public:
  vect::const_view<vec3r> r;
  sbox::pbc<real> const *simul_box;
  nl::data *nl;
  vect::const_view<nat_cont> all_contacts;
  vect::vector<nat_cont> *contacts;
  real cutoff;

public:
  void operator()() const;
  void mark_as_taken(nl::data *nl_) const;
};
} // namespace cg::nat_cont