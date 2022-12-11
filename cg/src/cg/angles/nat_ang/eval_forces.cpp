#include <cg/angles/nat_ang/eval_forces.h>

namespace cg::nat_ang {
void eval_forces::iter(int idx) const {
  auto angle = angles[idx];
  int i1 = angle.i1(), i2 = angle.i2(), i3 = angle.i3();
  auto nat_theta = angle.nat_theta();
  auto r1 = r[i1], r2 = r[i2], r3 = r[i3];

  auto r12 = r2 - r1, r23 = r3 - r2;

  auto x12_23 = cross(r12, r23);
  auto r12_rn = norm_inv(r12), r23_rn = norm_inv(r23);

  vec3r dtheta_dr1 = unit(cross(r12, x12_23)) * r12_rn;
  vec3r dtheta_dr3 = unit(cross(r23, x12_23)) * r23_rn;
  vec3r dtheta_dr2 = -dtheta_dr1 - dtheta_dr3;

  auto cos_theta = -dot(r12, r23) * r12_rn * r23_rn;
  cos_theta = clamp(cos_theta, (real)-1.0, (real)1.0);
  if (cos_theta == (real)0.0)
    dtheta_dr1 = dtheta_dr2 = dtheta_dr3 = vec3r::Zero();
  auto theta = acos(cos_theta);

  auto dtheta = theta - nat_theta;
  *V += CBA * dtheta * dtheta;

  auto dV_dtheta = (real)2.0 * CBA * dtheta;
  F[i1] -= dV_dtheta * dtheta_dr1;
  F[i2] -= dV_dtheta * dtheta_dr2;
  F[i3] -= dV_dtheta * dtheta_dr3;
}

int eval_forces::size() const {
  return angles.size();
}

} // namespace cg::nat_ang