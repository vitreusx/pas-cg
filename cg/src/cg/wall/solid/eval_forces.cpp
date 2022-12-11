#include <cg/wall/solid/eval_forces.h>

namespace cg::wall::solid {
void eval_forces::iter(int res_idx) const {
  auto res_r = r[res_idx];
  for (int wall_idx = 0; wall_idx < walls.size(); ++wall_idx) {
    auto const &wall = walls[wall_idx];
    auto sdist = wall.plane.signed_dist(res_r), dist = abs(sdist);
    if (dist >= min_dist)
      continue;

    auto dV_dr = depth * pow(dist, -10.0) * sign(sdist);
    F[res_idx] += dV_dr * wall.plane.normal();
    wall_F[wall_idx] -= dV_dr * wall.plane.normal();
  }
}

int eval_forces::size() const {
  return r.size();
}
} // namespace cg::wall::solid