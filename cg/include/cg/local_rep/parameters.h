#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::local_rep {
struct parameters {
  bool enabled;
  quantity r_excl, depth;

  void load(ioxx::xyaml::node const &n);
};
} // namespace cg::local_rep