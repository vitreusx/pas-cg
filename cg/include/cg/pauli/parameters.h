#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/xyaml.h>

namespace cg::pauli {
struct parameters {
  bool enabled;
  quantity r_excl, depth;
  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg::pauli