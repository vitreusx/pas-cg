#pragma once
#include "nat_cont.h"
#include <cg/nl/data.h>
#include <cg/types/box.h>

namespace cg::nat_cont {
class update_contacts {
public:
  nitro::const_view<vec3r> r;
  box<real> const *box;
  nl::data const *nl;
  nitro::const_view<nat_cont> all_contacts;
  nitro::vector<nat_cont> *contacts;

public:
  void operator()() const;
};
} // namespace cg::nat_cont