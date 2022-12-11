#include <cg/wall/lj/sift_free.h>

namespace cg::wall::lj {
void sift_free::iter(int res_idx) const {
  if (is_connected[res_idx])
    return;

  auto res_r = r[res_idx];
  for (int wall_idx = 0; wall_idx < walls.size(); ++wall_idx) {
    auto const &wall = walls[wall_idx];
    auto sdist = wall.plane.signed_dist(res_r), dist = abs(sdist);
    if (dist > min_dist)
      continue;

    auto [V_, dV_dr] = force(1.0 / dist);
    *V += V_;
    auto F_ = dV_dr * wall.plane.normal() * sign(sdist);
    F[res_idx] -= F_;
    wall_F[wall_idx] += F_;

    if (wall.num_conns < wall.limit) {
#pragma omp critical
      candidates->emplace_back(res_idx, wall_idx, dist);
    }
  }
}

int sift_free::size() const {
  return r.size();
}
} // namespace cg::wall::lj