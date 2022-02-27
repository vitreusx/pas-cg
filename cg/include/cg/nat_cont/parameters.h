#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::nat_cont {
struct parameters {
  bool enabled;
  quantity lj_depth, active_thr, ss_H1, ss_eqdist;
  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::nat_cont