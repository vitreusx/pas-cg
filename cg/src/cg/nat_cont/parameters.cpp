#include <cg/nat_cont/parameters.h>

namespace cg::nat_cont {
void parameters::unfolding_study_t::load(const ioxx::xyaml::node &n) {
  n["stop when all are formed/broken"] >> early_stopping;
  n["measure median times"] >> measure_times;
}

void parameters::load(ioxx::xyaml::node const &p) {
  p["enabled"] >> enabled;
  p["lj depth"] >> lj_depth;
  p["disulfide bond force"] >> ss_force;
  p["(un)folding study"] >> unfolding_study;
}
} // namespace cg::nat_cont