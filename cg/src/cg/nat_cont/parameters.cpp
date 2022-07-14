#include <cg/nat_cont/parameters.h>


namespace cg::nat_cont {

void parameters::load(ioxx::xyaml::node const &p) {
  p["enabled"] >> enabled;
  p["lj depth"] >> lj_depth;
  p["active threshold"] >> active_thr;
  p["disulfide bond force"] >> ss_force;
}
} // namespace cg::nat_cont