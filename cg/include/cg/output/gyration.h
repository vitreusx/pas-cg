#pragma once
#include "hook.h"
#include <cg/types/amp.h>

namespace cg::out {
class report_gyration_stuff : public hook {
public:
  nitro::const_view<vec3r> r;

public:
  void report_to(report_data &report) const override;
};
} // namespace cg::out