#include "heur/ang/eval_forces.h"
using namespace cg::heur_ang;

void eval_forces::operator()() const {
  for (int idx = 0; idx < angles.size(); ++idx) {
    iter(angles[idx]);
  }
}

template <typename E>
void eval_forces::iter(heur_ang_expr<E> const &angle) const {
  auto i1 = angle.i1(), i2 = angle.i2(), i3 = angle.i3();
  auto type_val = (uint8_t)angle.type();

  auto r1 = r[i1], r2 = r[i2], r3 = r[i3];
  auto r12 = r2 - r1, r23 = r3 - r2;

  auto x12_23 = cross(r12, r23);
  auto r12_rn = norm_inv(r12), r23_rn = norm_inv(r23);

  auto dtheta_dr1 = unit(cross(r12, x12_23)) * r12_rn;
  auto dtheta_dr3 = unit(cross(r23, x12_23)) * r23_rn;
  auto dtheta_dr2 = -dtheta_dr1 - dtheta_dr3;

  auto cos_theta = -dot(r12, r23) * r12_rn * r23_rn;
  auto theta = acos(cos_theta);

  real angle_V = 0.0f, dV_dtheta = 0.0f;
  for (int d = POLY_DEG; d >= 0; --d) {
    auto coeff = poly_coeffs[d][type_val];
    if (d > 0)
      dV_dtheta = (real)d * coeff + theta * dV_dtheta;
    angle_V = coeff + theta * angle_V;
  }

  *V += angle_V;

  F[i1] -= dV_dtheta * dtheta_dr1;
  F[i2] -= dV_dtheta * dtheta_dr2;
  F[i3] -= dV_dtheta * dtheta_dr3;
}

void eval_forces::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < angles.size(); ++idx) {
    iter(angles[idx]);
  }
}