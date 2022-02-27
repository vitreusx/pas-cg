#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::nl {
struct parameters {
  enum algorithm_t { LEGACY, CELL };
  algorithm_t algorithm;
  double pad_factor;
  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::nl