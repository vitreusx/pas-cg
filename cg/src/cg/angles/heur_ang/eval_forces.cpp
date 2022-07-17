#include <cg/angles/heur_ang/eval_forces.h>
namespace cg::heur_ang {

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
  cos_theta = clamp(cos_theta, (real)-1.0, (real)1.0);
  auto theta = acos(cos_theta);

  auto theta2 = theta * theta, theta3 = theta2 * theta, theta4 = theta3 * theta,
       theta5 = theta4 * theta, theta6 = theta3 * theta3;
  auto angle_V =
      poly_coeffs[0][type_val] + poly_coeffs[1][type_val] * theta +
      poly_coeffs[2][type_val] * theta2 + poly_coeffs[3][type_val] * theta3 +
      poly_coeffs[4][type_val] * theta4 + poly_coeffs[5][type_val] * theta5 +
      poly_coeffs[6][type_val] * theta6;
  auto dV_dtheta = poly_coeffs[1][type_val] +
                   (2.0 * poly_coeffs[2][type_val]) * theta +
                   (3.0 * poly_coeffs[3][type_val]) * theta2 +
                   (4.0 * poly_coeffs[4][type_val]) * theta3 +
                   (5.0 * poly_coeffs[5][type_val]) * theta4 +
                   (6.0 * poly_coeffs[6][type_val]) * theta5;

  *V += angle_V;

  //  dV_dtheta = clamp(dV_dtheta, (real)-1.0e3, (real)1.0e3);
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
} // namespace cg::heur_ang