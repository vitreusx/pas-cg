#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::tether {
struct parameters {
  bool enabled;
  quantity H1, H2, def_length;
  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::tether