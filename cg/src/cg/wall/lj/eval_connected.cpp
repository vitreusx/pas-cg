#include <cg/wall/lj/eval_connected.h>

namespace cg::wall::lj {
void eval_connected::iter(int conn_idx) const {
  decltype(auto) node = conns->at(conn_idx);
  if (node.is_vacant())
    return;

  auto &conn = node.item();
  auto &wall = walls[conn.wall_idx()];

  auto res_r = r[conn.res_idx()];
  auto bead_r = wall.plane.origin() + conn.bead_offset();
  auto r12 = bead_r - res_r;
  auto r12_n = norm(r12), r12_rn = (real)1.0 / r12_n;
  auto r12_u = r12 * r12_rn;

  auto [V_, dV_dr] = force(r12_rn);
  *V += conn.saturation() * V_;
  auto F_ = conn.saturation() * dV_dr * r12_u;
  F[conn.res_idx()] -= F_;
  wall_F[conn.wall_idx()] += F_;

  if (r12_n < breaking_dist) {
    conn.saturation() = min(conn.saturation() + saturation_diff, (real)1.0);
  } else {
    conn.saturation() = max(conn.saturation() - saturation_diff, (real)0.0);
    if (conn.saturation() == (real)0.0) {
#pragma omp critical
      removed->push_back(conn_idx);
    }
  }
}

int eval_connected::size() const {
  return conns->size();
}

} // namespace cg::wall::lj