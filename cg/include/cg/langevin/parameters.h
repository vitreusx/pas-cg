#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>
#include <variant>

namespace cg::lang {
struct parameters {
  bool enabled;
  quantity gamma, dt;
  using quantity_range = std::pair<quantity, quantity>;
  std::variant<quantity, quantity_range> temperature;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::lang