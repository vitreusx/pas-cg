#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/xyaml.h>

namespace cg::cnd {
struct parameters {
  bool enabled;
  quantity CDA, CDB;

  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg::cnd