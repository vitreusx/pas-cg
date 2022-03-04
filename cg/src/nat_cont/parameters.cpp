#include "nat_cont/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::nat_cont;

void parameters::load(ioxx::xyaml::node const &p) {
  p["enabled"] >> enabled;
  p["lj depth"] >> lj_depth;
  p["active threshold"] >> active_thr;
  p["disulfide bond force"] >> disulfide_force;
}