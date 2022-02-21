#include "features/simp_nat_dih/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::snd;

void parameters::connect(ioxx::xyaml_proxy &p) {
  enabled = p["enabled"].as<bool>();
  CDH = p["CDH"].as<quantity>();
}