#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/xyaml.h>

namespace cg::nat_cont {
struct parameters {
  bool enabled;
  quantity lj_depth, active_thr, ss_H1, ss_eqdist;
  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg::nat_cont