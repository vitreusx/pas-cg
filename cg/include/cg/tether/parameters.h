#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::tether {
struct parameters {
  bool enabled;
  quantity H1, H2, def_length;
  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::tether