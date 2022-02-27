#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::const_dh {
struct parameters {
  bool enabled;
  quantity screening_dist, permittivity;
  void link(ioxx::xyaml::proxy &proxy);
};
} // namespace cg::const_dh