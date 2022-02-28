#pragma once
#include "nat_cont.h"
#include "parameters.h"
#include <cg/output/hook.h>

namespace cg::nat_cont {
class report_stuff : public out::hook {
public:
  nitro::const_view<nat_cont> all_nat_conts;
  nitro::const_view<vec3r> r;
  parameters const *params;

public:
  void report_to(out::report_state &report) const override;
};
} // namespace cg::nat_cont