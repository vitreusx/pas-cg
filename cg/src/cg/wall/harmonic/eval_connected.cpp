#include <cg/base_forces/harmonic.h>
#include <cg/wall/harmonic/eval_connected.h>

namespace cg::wall::harmonic {
void eval_connected::iter(int idx) const {
  auto conn = conns[idx];
  auto res_r = r[conn.res_idx()];
  auto bead_r = conn.bead_offset() + walls[conn.wall_idx()].plane.origin();
  auto r12 = bead_r - res_r;
  auto r12_n = norm(r12), r12_rn = (real)1.0 / r12_n;
  auto r12_u = r12 * r12_rn;

  auto [V_, dV_dr] = cg::harmonic(HH1, (real)0.0, conn.nat_dist())(r12_n);
  *V += V_;
  F[conn.res_idx()] -= dV_dr * r12_u;
  wall_F[conn.wall_idx()] += dV_dr * r12_u;
}

int eval_connected::size() const {
  return conns.size();
}
} // namespace cg::wall::harmonic