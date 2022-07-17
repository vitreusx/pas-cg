#include <cg/output/print_raw_data.h>

namespace cg::out {
void print_raw_data::operator()() const {
  *data_file << std::setprecision(17);
  *data_file << st->traj_idx << " " << st->step_idx << " " << st->dyn.V << " "
             << st->num_res << '\n';
  for (int idx = 0; idx < st->num_res; ++idx) {
    auto r = st->r[idx], F = st->dyn.F[idx];
    *data_file << r.x() << " " << r.y() << " " << r.z() << " " << F.x() << " "
               << F.y() << " " << F.z() << '\n';
  }
}
} // namespace cg::out