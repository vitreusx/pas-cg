#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/xyaml.h>

namespace cg::chir {
struct parameters {
  bool enabled;
  quantity e_chi;
  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg::chir