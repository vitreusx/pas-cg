#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::chir {
struct parameters {
  bool enabled;
  quantity e_chi;
  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::chir