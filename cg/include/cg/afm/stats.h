#pragma once
#include "compiled.h"
#include "vel/eval_forces.h"
#include <cg/output/hook.h>

namespace cg::afm {
class report_stats : public out::hook {
public:
  compiled_tips const *tips;
  nitro::const_view<vec3r> r, v;
  vel::eval_forces const *eval_vel_forces;
  nitro::const_view<int> chain_first, chain_last;

public:
  void report_to(out::report_state &report) const override;
};
} // namespace cg::afm