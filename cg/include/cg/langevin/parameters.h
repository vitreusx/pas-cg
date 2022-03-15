#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::lang {
enum class lang_type { LEGACY, NORMAL };

struct parameters {
  bool enabled;
  quantity gamma, dt;
  using quantity_range = std::pair<quantity, quantity>;
  std::variant<quantity, quantity_range> temperature;
  lang_type type;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::lang