#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::rel_dh {
struct parameters {
  bool enabled;
  quantity screening_dist, perm_factor;
  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::rel_dh