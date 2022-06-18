#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::nl {
struct parameters {
  enum algorithm_t {
    LEGACY,
    CELL
  };
  algorithm_t algorithm;
  quantity pad;
  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::nl