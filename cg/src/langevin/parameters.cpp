#include "langevin/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::lang;

void parameters::connect(ioxx::xyaml_proxy &p) {
  enabled = p["enabled"].as<bool>();
  gamma = p["gamma factor"].as<quantity>();
  temperature = p["temperature"].as<quantity>();
  dt = p["dt"].as<quantity>();
}