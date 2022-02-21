#include "features/chir/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::chir;

void parameters::connect(ioxx::xyaml_proxy &p) {
  enabled = p["enabled"].as<bool>();
  e_chi = p["e_chi"].as<quantity>();
}