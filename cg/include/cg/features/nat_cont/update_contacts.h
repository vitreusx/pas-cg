#pragma once
#include "nat_cont.h"
#include <cg/features/nl/data.h>
#include <cg/types/box.h>

namespace cg::nat_cont {
class update_contacts {
public:
  nitro::vector<vec3r> const *r;
  box<real> const *box;
  nl::data const *nl;
  nitro::vector<nat_cont> const *all_contacts;
  nitro::vector<nat_cont> *contacts;

public:
  void operator()() const;
};
} // namespace cg::nat_cont