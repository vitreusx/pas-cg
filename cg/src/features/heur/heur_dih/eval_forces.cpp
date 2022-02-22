#include "features/heur/heur_dih/eval_forces.h"
using namespace cg::heur_dih;

void eval_forces::operator()() const {
  for (int idx = 0; idx < dihedrals->size(); ++idx) {
    iter(dihedrals->at(idx));
  }
}

template <typename E>
void eval_forces::iter(heur_dih_expr<E> const &dihedral) const {
  auto i1 = dihedral.i1(), i2 = dihedral.i2(), i3 = dihedral.i3(),
       i4 = dihedral.i4();
  auto type_val = (int8_t)dihedral.type();

  auto r1 = r->at(i1), r2 = r->at(i2), r3 = r->at(i3), r4 = r->at(i4);
  auto r12 = r2 - r1, r23 = r3 - r2, r34 = r4 - r3;
  auto x12_23 = cross(r12, r23), x23_34 = cross(r23, r34);

  auto x12_23_rn = norm_inv(x12_23), x23_34_rn = norm_inv(x23_34);
  auto x12_23_u = x12_23 * x12_23_rn, x23_34_u = x23_34 * x23_34_rn;

  auto cos_phi = dot(x12_23_u, x23_34_u);
  auto phi = cos(cos_phi);
  if (dot(x12_23, r34) < 0.0f)
    phi = -phi;
  auto sin_phi = sin(phi);

  auto sin2_phi = sin_phi * sin_phi, cos2_phi = cos_phi * cos_phi,
       sin_phi_cos_phi = sin_phi * cos_phi;

  *V += coeffs[0][type_val] + coeffs[1][type_val] * sin_phi +
        coeffs[2][type_val] * cos_phi + coeffs[3][type_val] * sin2_phi +
        coeffs[4][type_val] * cos2_phi + coeffs[5][type_val] * sin_phi_cos_phi;

  auto dV_dphi =
      coeffs[1][type_val] * cos_phi - coeffs[2][type_val] * sin_phi +
      2.0f * (coeffs[3][type_val] + coeffs[4][type_val]) * sin_phi_cos_phi +
      coeffs[5][type_val] * (cos2_phi - sin2_phi);

  auto r23_n = norm(r23);
  auto dphi_dr1 = -x12_23_u * r23_n * x12_23_rn;
  auto dphi_dr4 = x23_34_u * r23_n * x23_34_rn;
  auto df =
      (-dphi_dr1 * dot(r12, r23) + dphi_dr4 * dot(r23, r34)) / (r23_n * r23_n);
  auto dphi_dr2 = -dphi_dr1 + df;
  auto dphi_dr3 = -dphi_dr4 - df;

  F->at(i1) -= dV_dphi * dphi_dr1;
  F->at(i2) -= dV_dphi * dphi_dr2;
  F->at(i3) -= dV_dphi * dphi_dr3;
  F->at(i4) -= dV_dphi * dphi_dr4;
}

void eval_forces::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < dihedrals->size(); ++idx) {
    iter(dihedrals->at(idx));
  }
}