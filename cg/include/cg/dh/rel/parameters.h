#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/xyaml.h>

namespace cg::rel_dh {
struct parameters {
  bool enabled;
  quantity screening_dist, perm_factor;
  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg::rel_dh