#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::nat_ang {
struct parameters {
  bool enabled;
  quantity CBA;
  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::nat_ang