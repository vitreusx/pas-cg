#pragma once
#include <cg/base_forces/disulfide.h>
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::nat_cont {
struct parameters {
  bool enabled;
  quantity lj_depth, active_thr;
  std::optional<disulfide_force> disulfide_force;
  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::nat_cont