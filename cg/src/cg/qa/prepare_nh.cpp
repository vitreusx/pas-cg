#include <cg/qa/prepare_nh.h>

namespace cg::qa {
void prepare_nh::iter(int idx) const {
  int iprev = prev[idx], icur = idx, inext = next[idx];
  if (iprev < 0 || inext < 0)
    return;

  auto rprev = r[iprev], rcur = r[icur], rnext = r[inext];
  auto v1 = rcur - rprev, v2 = rnext - rcur;
  n[icur] = unit(v2 - v1);
  h[icur] = unit(cross(v2, v1));
}

int prepare_nh::size() const {
  return num_particles;
}
} // namespace cg::qa