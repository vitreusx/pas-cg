#pragma once
#include "hook.h"
#include <cg/types/amp.h>

namespace cg::out {
class compute_rmsd : public hook {
public:
  nitro::const_view<vec3r> orig_r, cur_r;

public:
  void report_to(report_data &report) const override;
};
} // namespace cg::out