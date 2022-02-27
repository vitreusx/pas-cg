#pragma once
#include <cg/types/vec3.h>
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::fafm {
struct parameters {
  bool enabled;

  struct afm_tip_t {
    int pulled_idx;
    vec3<quantity> F;
  };
  std::vector<afm_tip_t> afm_tips;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::fafm