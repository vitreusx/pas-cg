#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::pauli {
struct parameters {
  bool enabled;
  quantity r_excl, depth;
  void link(ioxx::xyaml::proxy &proxy);
};
} // namespace cg::pauli