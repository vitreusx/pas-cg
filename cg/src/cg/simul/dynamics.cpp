#include <cg/simul/dynamics.h>
#include <cg/simul/thread.h>
namespace cg::simul {

dynamics::dynamics(int num_residues) {
  F = vect::vector<vec3r>(num_residues);
}

void dynamics::reset() {
  V = 0.0;
  for (int idx = 0; idx < F.size(); ++idx)
    F[idx] = vec3r::Zero();
}

void dynamics::omp_reset() {
#pragma omp master
  V = 0.0;

#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < F.size(); ++idx)
    F[idx] = vec3r::Zero();
}

void dynamics::reduce(dynamics &target) {
  target.V += V;

  for (int idx = 0; idx < F.size(); ++idx)
    target.F[idx] += F[idx];
}

void dynamics::omp_reduce(dynamics &target) {
#pragma omp atomic
  target.V += V;

  for (int idx = 0; idx < F.size(); ++idx)
    target.F[idx].omp_atomic_add(F[idx]);
}

void dynamics::omp_reduce_v2(dynamics &target, thread const &thr) {
  auto num_thr = thr.team->num_threads;
  auto N = F.size();

  for (int tid = 0; tid < num_thr; ++tid) {
    if (thr.tid == tid)
      target.V += V;

    auto sector = (tid + thr.tid) % num_thr;
    auto start = (N * sector) / num_thr, end = (N * (sector + 1)) / num_thr;
    for (int idx = start; idx < end; ++idx)
      target.F[idx] += F[idx];

#pragma omp barrier
  }
}

void dynamics::omp_reduce_v3(dynamics &target, const thread_team &team) {
#pragma omp atomic
  target.V += V;

#pragma omp barrier

#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < F.size(); ++idx) {
    vec3r f;
    for (auto const &F_ : team.forces)
      f += F_[idx];
    target.F[idx] += f;
  }
}

} // namespace cg::simul