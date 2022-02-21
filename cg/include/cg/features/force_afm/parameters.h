#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/xyaml.h>
#include <cg/types/vec3.h>

namespace cg::fafm {
struct parameters {
  bool enabled;

  struct afm_tip_t {
    int pulled_idx;
    vec3<quantity> F;
  };
  std::vector<afm_tip_t> afm_tips;

  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg::fafm