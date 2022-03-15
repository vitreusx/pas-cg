#pragma once
#include "hook.h"
#include <cg/input/model.h>
#include <cg/types/amp.h>

namespace cg::out {
class add_structure : public hook {
public:
  input::model const *model;
  input::model::res_map_t const *res_map;
  nitro::const_view<vec3r> r;

  void report_to(report_data &report) const override;
};
} // namespace cg::out