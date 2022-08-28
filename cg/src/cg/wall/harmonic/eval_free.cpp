#include <cg/wall/harmonic/eval_free.h>

namespace cg::wall::harmonic {
void eval_free::operator()() const {
  for (int idx = 0; idx < r.size(); ++idx)
    iter(idx);
}

void eval_free::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < r.size(); ++idx)
    iter(idx);
}

void eval_free::iter(int idx) const {
  if (!is_connected[idx]) {
    auto res_r = r[idx];
    for (int wall_idx = 0; wall_idx < walls.size(); ++wall_idx) {
      auto const &wall = walls[wall_idx];
      auto sdist = wall.plane.signed_dist(res_r), dist = abs(sdist);
      if (dist >= min_dist)
        continue;

      auto dV_dr = depth * pow(dist, -10.0) * sign(sdist);
      auto force = dV_dr * wall.plane.normal();
      F[idx] -= force;
      wall_F[wall_idx] += force;
    }
  }
}

void eval_free::for_slice(int from, int to) const {
  for (int idx = from; idx < to; ++idx)
    iter(idx);
}

int eval_free::total_size() const {
  return r.size();
}


} // namespace cg::wall::harmonic