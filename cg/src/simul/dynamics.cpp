#include "simul/dynamics.h"
namespace cg::simul {

dynamics::dynamics(int num_residues) { F = nitro::vector<vec3r>(num_residues); }

void dynamics::reset() {
  V = 0.0;
  for (int idx = 0; idx < F.size(); ++idx)
    F[idx] = vec3r::Zero();
}

void dynamics::omp_reset() {
#pragma omp atomic write
  V = 0.0;

#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < F.size(); ++idx)
    F[idx] = vec3r::Zero();
}

void dynamics::omp_reduce(dynamics &target) {
#pragma omp atomic
  target.V += V;

  for (int idx = 0; idx < F.size(); ++idx)
    target.F[idx].omp_atomic_add(F[idx]);
}

void dynamics::sanity_check() {
  if (isnan(V))
    throw;

  for (int idx = 0; idx < F.size(); ++idx) {
    if (isnan(F[idx].x()))
      throw;
    if (isnan(F[idx].y()))
      throw;
    if (isnan(F[idx].z()))
      throw;
  }
}
} // namespace cg::simul