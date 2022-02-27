#include "dh/rel/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::rel_dh;

void parameters::load(ioxx::xyaml::node const &p) {
  enabled = p["enabled"].as<bool>();
  screening_dist = p["screening distance"].as<quantity>();
  perm_factor = p["factor"].as<quantity>();
}