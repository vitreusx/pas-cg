#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::pauli {
struct parameters {
  bool enabled;
  quantity r_excl, depth;
  void link(ioxx::xyaml::proxy &proxy);
};
} // namespace cg::pauli