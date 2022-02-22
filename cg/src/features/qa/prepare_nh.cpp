#include "features/qa/prepare_nh.h"
using namespace cg::qa;

void prepare_nh::operator()() const {
  for (int idx = 0; idx < num_particles; ++idx) {
    iter(idx);
  }
}

void prepare_nh::iter(int idx) const {
  auto iprev = prev->at(idx), icur = idx, inext = next->at(idx);
  if (iprev < 0 || inext < 0)
    return;

  auto rprev = r->at(iprev), rcur = r->at(icur), rnext = r->at(inext);
  auto v1 = rcur - rprev, v2 = rnext - rcur;
  n->at(icur) = unit(v2 - v1);
  h->at(icur) = unit(cross(v2, v1));
}

void prepare_nh::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < num_particles; ++idx) {
    iter(idx);
  }
}