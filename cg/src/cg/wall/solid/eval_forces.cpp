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

void eval_forces::operator()() const {
  for (int idx = 0; idx < r.size(); ++idx)
    iter(idx);
}

void eval_forces::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < r.size(); ++idx)
    iter(idx);
}

void eval_forces::for_slice(int from, int to) const {
  for (int idx = from; idx < to; ++idx)
    iter(idx);
}

int eval_forces::total_size() const {
  return r.size();
}


} // namespace cg::wall::solid