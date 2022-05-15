#pragma once
#include <cg/base_forces/disulfide.h>
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::nat_cont {
struct parameters {
  bool enabled;
  quantity lj_depth, active_thr;
  std::optional<disulfide_force> ss_force;
  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::nat_cont