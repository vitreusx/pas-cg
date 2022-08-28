#include <cg/angles/heur_dih/rel.h>
namespace cg::heur_dih {

void eval_forces::operator()() const {
  for (int idx = 0; idx < dihedrals.size(); ++idx) {
    iter(dihedrals[idx]);
  }
}

template <typename E>
void eval_forces::iter(heur_dih_expr<E> const &dihedral) const {
  auto i1 = dihedral.i1(), i2 = dihedral.i2(), i3 = dihedral.i3(),
       i4 = dihedral.i4();
  auto type_val = (int8_t)dihedral.type();

  auto r1 = r[i1], r2 = r[i2], r3 = r[i3], r4 = r[i4];
  auto r12 = r2 - r1, r23 = r3 - r2, r34 = r4 - r3;
  auto x23_12 = cross(r23, r12), x34_23 = cross(r34, r23);

  auto x23_12_rn = norm_inv(x23_12), x34_23_rn = norm_inv(x34_23);
  auto x23_12_u = x23_12 * x23_12_rn, x34_23_u = x34_23 * x34_23_rn;

  auto cos_phi = dot(x23_12_u, x34_23_u);
  auto phi = acos(cos_phi);
  if (-dot(r12, x34_23_u) < 0.0f)
    phi = -phi;
  auto sin_phi = sin(phi);

  auto sin2_phi = sin_phi * sin_phi, cos2_phi = cos_phi * cos_phi,
       sin_phi_cos_phi = sin_phi * cos_phi;

  auto V_ = coeffs.const_[type_val] + coeffs.sin[type_val] * sin_phi +
            coeffs.cos[type_val] * cos_phi + coeffs.sin2[type_val] * sin2_phi +
            coeffs.cos2[type_val] * cos2_phi +
            coeffs.sin_cos[type_val] * sin_phi_cos_phi;
  *V += V_;

  auto dV_dphi =
      coeffs.sin[type_val] * cos_phi - coeffs.cos[type_val] * sin_phi +
      2.0f * (coeffs.sin2[type_val] - coeffs.cos2[type_val]) * sin_phi_cos_phi +
      coeffs.sin_cos[type_val] * (cos2_phi - sin2_phi);

  auto r23_n = norm(r23);
  auto dphi_dr1 = x23_12_u * r23_n * x23_12_rn;
  auto dphi_dr4 = -x34_23_u * r23_n * x34_23_rn;
  auto df = (dphi_dr1 * (-dot(r12, r23)) - dphi_dr4 * (-dot(r34, r23))) /
            (r23_n * r23_n);
  auto dphi_dr2 = -dphi_dr1 + df;
  auto dphi_dr3 = -dphi_dr4 - df;

  F[i1] -= dV_dphi * dphi_dr1;
  F[i2] -= dV_dphi * dphi_dr2;
  F[i3] -= dV_dphi * dphi_dr3;
  F[i4] -= dV_dphi * dphi_dr4;
}

void eval_forces::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < dihedrals.size(); ++idx) {
    iter(dihedrals[idx]);
  }
}

void eval_forces::for_slice(int from, int to) const {
  for (int idx = from; idx < to; ++idx)
    iter(dihedrals[idx]);
}

int eval_forces::total_size() const {
  return dihedrals.size();
}


} // namespace cg::heur_dih