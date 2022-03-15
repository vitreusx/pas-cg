#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::nl {
struct parameters {
  enum algorithm_t { LEGACY, CELL };
  algorithm_t algorithm;
  quantity pad;
  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::nl