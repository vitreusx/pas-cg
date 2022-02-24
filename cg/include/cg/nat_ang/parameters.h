#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/xyaml.h>

namespace cg::nat_ang {
struct parameters {
  bool enabled;
  quantity k;
  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg::nat_ang