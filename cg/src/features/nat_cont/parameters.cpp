#include "features/nat_cont/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::nat_cont;

void parameters::connect(ioxx::xyaml_proxy &p) {
  enabled = p["enabled"].as<bool>();
  lj_depth = p["lj depth"].as<quantity>();
  active_thr = p["active threshold"].as<double>();
  ss_H1 = p["disulfide bonds"]["H1"].as<quantity>();
  ss_eqdist = p["disulfide bonds"]["equilibrium dist"].as<quantity>();
}