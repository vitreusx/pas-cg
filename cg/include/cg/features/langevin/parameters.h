#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/xyaml.h>

namespace cg::lang {
struct parameters {
  bool enabled;
  quantity gamma, temperature, dt;

  void connect(ioxx::xyaml_proxy& proxy);
};
} // namespace cg::lang