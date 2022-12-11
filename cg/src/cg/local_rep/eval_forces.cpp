#include <cg/base_forces/shifted_lj.h>
#include <cg/local_rep/eval_forces.h>

namespace cg::local_rep {
void eval_forces::iter(int idx) const {
  auto pair = pairs[idx];
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

int eval_forces::size() const {
  return pairs.size();
}

} // namespace cg::local_rep