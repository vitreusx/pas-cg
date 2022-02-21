#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/xyaml.h>

namespace cg::tether {
struct parameters {
  bool enabled;
  quantity H1, H2, def_length;
  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg::tether