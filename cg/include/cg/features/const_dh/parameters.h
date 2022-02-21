#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/xyaml.h>

namespace cg::const_dh {
struct parameters {
  bool enabled;
  quantity screening_dist, permittivity;
  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg::const_dh