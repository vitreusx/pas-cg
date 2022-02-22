#include "features/comp_nat_dih/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::cnd;

void parameters::connect(ioxx::xyaml_proxy &p) {
  enabled = p["enabled"].as<bool>();
  CDA = p["CDA"].as<quantity>();
  CDB = p["CDB"].as<quantity>();
}