#include <cg/base_forces/shifted_lj.h>
#include <cg/local_rep/eval_forces.h>

namespace cg::local_rep {
void eval_forces::operator()() const {
  for (int idx = 0; idx < pairs.size(); ++idx) {
    iter(pairs[idx]);
  }
}

template <typename E> void eval_forces::iter(pair_expr<E> const &pair) const {
  auto i1 = pair.i1(), i2 = pair.i2();

  auto r1 = r[i1], r2 = r[i2];
  auto r12 = r2 - r1;
  auto r12_rn = norm_inv(r12);

  if (1.0f > r12_rn * cutoff)
    return;

  auto r12_u = r12 * r12_rn;
  auto [V_, dV_dr] = shifted_lj(depth, cutoff)(r12_rn);
  dV_dr = clamp(dV_dr, (real)-1e3, (real)1e3);

  *V += V_;
  F[i1] += r12_u * dV_dr;
  F[i2] -= r12_u * dV_dr;
}

void eval_forces::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < pairs.size(); ++idx) {
    iter(pairs[idx]);
  }
}
} // namespace cg::local_rep