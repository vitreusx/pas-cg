#include <algorithm>
#include <cg/wall/lj/process_candidates.h>

namespace cg::wall::lj {
void process_candidates::operator()() const {
  for (int idx = 0; idx < removed->size(); ++idx) {
    auto conn_idx = removed->at(idx);
    auto conn = conns->at(conn_idx).item();
    is_connected[conn.res_idx()] = false;
    --walls[conn.wall_idx()].num_conns;
    conns->remove(conn_idx);
  }
  removed->clear();

  std::sort(candidates->begin(), candidates->end(),
            [](auto const &x, auto const &y) -> bool {
              return x.dist() < y.dist();
            });

  for (int idx = 0; idx < candidates->size(); ++idx) {
    auto const &cand = candidates->at(idx);
    auto &wall = walls[cand.wall_idx()];
    if (is_connected[cand.res_idx()] || wall.num_conns >= wall.limit)
      continue;

    auto res_r = r[cand.res_idx()];
    auto bead_offset = wall.plane.projection(res_r) - wall.plane.origin();
    auto saturation = 0.0f;
    conns->emplace(cand.res_idx(), cand.wall_idx(), bead_offset, saturation);
    ++wall.num_conns;
  }
}
} // namespace cg::wall::lj