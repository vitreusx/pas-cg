#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::cnd {
struct parameters {
  bool enabled;
  quantity CDA, CDB;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::cnd