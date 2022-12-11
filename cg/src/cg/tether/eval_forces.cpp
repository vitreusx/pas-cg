#include <cg/base_forces/harmonic.h>
#include <cg/tether/eval_forces.h>

namespace cg::tether {
void eval_forces::iter(int idx) const {
  auto tether = tethers[idx];
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

int eval_forces::size() const {
  return tethers.size();
}
} // namespace cg::tether