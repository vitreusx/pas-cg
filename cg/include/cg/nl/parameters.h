#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::nl {
struct parameters {
  quantity pad, cutoff;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::nl