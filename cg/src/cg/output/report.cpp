#include <cg/output/report.h>
#include <limits>

namespace cg::out {

void report::traj_init(int traj_idx) {
  stats_last = std::numeric_limits<real>::lowest();
  struct_last = std::numeric_limits<real>::lowest();

  auto &traj_out_div = out_root.add<ioxx::sl4::div>();
  traj_out_div.add<ioxx::sl4::comment>("trajectory ", traj_idx);
  traj_out_div.named_add<ioxx::sl4::table>("scalars");

  auto &traj_map_div = map_root.add<ioxx::sl4::div>();
  traj_map_div.add<ioxx::sl4::comment>("trajectory ", traj_idx);
}

} // namespace cg::out