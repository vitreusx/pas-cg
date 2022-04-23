#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::const_dh {
struct parameters {
  bool enabled;
  quantity screening_dist, permittivity;
  void link(ioxx::xyaml::proxy &proxy);
};
} // namespace cg::const_dh