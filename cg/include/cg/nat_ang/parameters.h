#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::nat_ang {
struct parameters {
  bool enabled;
  quantity k;
  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::nat_ang