#include "pauli/eval_forces.h"
#include "base_forces/shifted_lj.h"
using namespace cg::pauli;

void eval_forces::operator()() const {
  for (int idx = 0; idx < pairs->size(); ++idx) {
    iter(pairs->at(idx));
  }
}

template <typename E> void eval_forces::iter(pair_expr<E> const &pair) const {
  auto i1 = pair.i1(), i2 = pair.i2();

  auto r1 = r[i1], r2 = r[i2];
  auto r12 = simul_box->r_uv(r1, r2);
  auto r12_rn = norm_inv(r12);

  if (1.0f < r12_rn * r_excl) {
    auto r12_u = r12 * r12_rn;
    auto [V_, dV_dr] = shifted_lj(depth, r_excl)(r12_rn);

    *V += V_;
    F[i1] += r12_u * dV_dr;
    F[i2] -= r12_u * dV_dr;
  }
}

void eval_forces::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < pairs->size(); ++idx) {
    iter(pairs->at(idx));
  }
}