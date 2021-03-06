#include <cg/wall/solid.h>

namespace cg::wall::solid {
void eval_forces::iter(int res_idx) const {
  auto res_r = r[res_idx];
  for (int wall_idx = 0; wall_idx < walls.size(); ++wall_idx) {
    auto const &wall = walls[wall_idx];
    auto sdist = wall.signed_dist(res_r), dist = abs(sdist);
    if (dist >= min_dist)
      continue;

    auto dV_dr = depth * pow(dist, -10.0) * sign(sdist);
    F[res_idx] += dV_dr * wall.normal();
    wall_F[wall_idx] -= dV_dr * wall.normal();
  }
}

void eval_forces::operator()() const {
  for (int idx = 0; idx < r.size(); ++idx)
    iter(idx);
}

void eval_forces::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < r.size(); ++idx)
    iter(idx);
}
} // namespace cg::wall::solid