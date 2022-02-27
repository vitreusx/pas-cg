#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::chir {
struct parameters {
  bool enabled;
  quantity e_chi;
  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::chir