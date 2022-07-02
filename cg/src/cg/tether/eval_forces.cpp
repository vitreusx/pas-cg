#include <cg/base_forces/harmonic.h>
#include <cg/tether/eval_forces.h>
namespace cg::tether {

void eval_forces::operator()() const {
  for (int idx = 0; idx < tethers.size(); ++idx) {
    iter(tethers[idx]);
  }
}

void eval_forces::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < tethers.size(); ++idx) {
    iter(tethers[idx]);
  }
}

template <typename E> void eval_forces::iter(pair_expr<E> const &tether) const {
  auto i1 = tether.i1(), i2 = tether.i2();
  auto nat_dist = tether.nat_dist();

  auto r1 = r[i1], r2 = r[i2];
  auto r12 = r2 - r1;

  auto r12_n = norm(r12);
  auto [V_, dV_dr] = harmonic(H1, H2, nat_dist)(r12_n);
  dV_dr = clamp(dV_dr, (real)-1e3, (real)1e3);

  auto r12_u = r12 / r12_n;
  *V += V_;
  F[i1] += dV_dr * r12_u;
  F[i2] -= dV_dr * r12_u;
}
} // namespace cg::tether