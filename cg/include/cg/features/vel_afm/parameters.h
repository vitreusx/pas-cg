#pragma once
#include <cg/types/vec3.h>
#include <cg/utils/quantity.h>
#include <ioxx/xyaml.h>

namespace cg::vafm {
struct parameters {
  bool enabled;
  quantity H1, H2;

  struct afm_tip_t {
    int pulled_idx;
    std::optional<vec3<quantity>> origin;
    vec3<quantity> v;
  };
  std::vector<afm_tip_t> afm_tips;

  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg::vafm