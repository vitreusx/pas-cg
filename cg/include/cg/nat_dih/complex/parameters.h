#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::cnd {
struct parameters {
  bool enabled;
  quantity CDA, CDB;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::cnd