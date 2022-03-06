#pragma once
#include "nat_cont.h"
#include "parameters.h"
#include <cg/output/hook.h>

namespace cg::nat_cont {
class report_stuff : public out::hook {
public:
  nitro::const_view<nat_cont> all_contacts;
  nitro::const_view<vec3r> r;
  nitro::const_view<int> chain_idx;
  parameters const *params;

  quantity cur_dist(nat_cont cont) const;
  bool is_active(nat_cont cont) const;

public:
  void report_to(out::report_state &report) const override;
};
} // namespace cg::nat_cont