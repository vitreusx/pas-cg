#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::lang {
struct parameters {
  bool enabled;
  quantity gamma, temperature, dt;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::lang