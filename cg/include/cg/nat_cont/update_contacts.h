#pragma once
#include "nat_cont.h"
#include <cg/nl/data.h>
#include <cg/types/box.h>

namespace cg::nat_cont {
class update_contacts {
public:
  vect::const_view<vec3r> r;
  box<real> const *simul_box;
  nl::data const *nl;
  vect::const_view<nat_cont> all_contacts;
  vect::vector<nat_cont> *contacts;

public:
  void operator()() const;
};
} // namespace cg::nat_cont