#include "nat_dih/complex/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::cnd;

void parameters::load(ioxx::xyaml::node const &p) {
  enabled = p["enabled"].as<bool>();
  CDA = p["CDA"].as<quantity>();
  CDB = p["CDB"].as<quantity>();
}