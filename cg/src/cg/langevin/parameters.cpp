#include <cg/langevin/parameters.h>

namespace cg::lang {

void parameters::load(ioxx::xyaml::node const &p) {
  enabled = p["enabled"].as<bool>();
  gamma = p["gamma factor"].as<quantity>();
  temperature = p["temperature"].as<quantity>();
  dt = p["dt"].as<quantity>();
}
} // namespace cg::lang