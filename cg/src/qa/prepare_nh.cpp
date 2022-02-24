#include "qa/prepare_nh.h"
using namespace cg::qa;

void prepare_nh::operator()() const {
  for (int idx = 0; idx < num_particles; ++idx) {
    iter(idx);
  }
}

void prepare_nh::iter(int idx) const {
  int iprev = prev[idx], icur = idx, inext = next[idx];
  if (iprev < 0 || inext < 0)
    return;

  auto rprev = r[iprev], rcur = r[icur], rnext = r[inext];
  auto v1 = rcur - rprev, v2 = rnext - rcur;
  n[icur] = unit(v2 - v1);
  h[icur] = unit(cross(v2, v1));
}

void prepare_nh::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < num_particles; ++idx) {
    iter(idx);
  }
}