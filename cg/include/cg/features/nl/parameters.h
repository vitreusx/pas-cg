#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/xyaml.h>

namespace cg::nl {
struct parameters {
  bool enabled;
  double pad_factor;
  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg::nl