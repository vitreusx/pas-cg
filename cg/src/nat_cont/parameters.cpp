#include "nat_cont/parameters.h"
#include "utils/ioxx_interop.h"
namespace cg::nat_cont {

void parameters::load(ioxx::xyaml::node const &p) {
  p["enabled"] >> enabled;
  p["lj depth"] >> lj_depth;
  p["active threshold"] >> active_thr;
  p["disulfide bond force"] >> ss_force;
}
} // namespace cg::nat_cont