#include "features/tether/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::tether;

void parameters::connect(ioxx::xyaml_proxy &p) {
  enabled = p["enabled"].as<bool>();
  H1 = p["H1"].as<quantity>();
  H2 = p["H2"].as<quantity>();
  def_length = p["default bond length"].as<quantity>();
}