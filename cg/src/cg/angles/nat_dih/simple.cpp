#include <cg/angles/nat_dih/simple.h>
namespace cg::snd {

void eval_forces::operator()() const {
  for (int idx = 0; idx < dihedrals.size(); ++idx) {
    iter(dihedrals[idx]);
  }
}

template <typename E>
void eval_forces::iter(nat_dih_expr<E> const &nat_dih) const {
  int i1 = nat_dih.i1(), i2 = nat_dih.i2(), i3 = nat_dih.i3(),
      i4 = nat_dih.i4();
  auto r1 = r[i1], r2 = r[i2], r3 = r[i3], r4 = r[i4];
  auto r12 = r2 - r1, r23 = r3 - r2, r34 = r4 - r3;
  auto x12_23 = cross(r12, r23), x23_34 = cross(r23, r34);

  auto x12_23_rn = norm_inv(x12_23), x23_34_rn = norm_inv(x23_34);
  auto x12_23_u = x12_23 * x12_23_rn, x23_34_u = x23_34 * x23_34_rn;

  auto cos_phi = dot(x12_23_u, x23_34_u);
  auto phi = acos(cos_phi);
  if (dot(x12_23, r34) < 0.0f)
    phi = -phi;

  auto diff = phi - nat_dih.nat_phi();
  *V += 0.5f * CDH * diff * diff;

  auto r23_n = norm(r23);
  auto dphi_dr1 = -x12_23_u * r23_n * x12_23_rn;
  auto dphi_dr4 = x23_34_u * r23_n * x23_34_rn;
  auto df =
      (-dphi_dr1 * dot(r12, r23) + dphi_dr4 * dot(r23, r34)) / (r23_n * r23_n);
  auto dphi_dr2 = -dphi_dr1 + df;
  auto dphi_dr3 = -dphi_dr4 - df;

  auto dV_dphi = CDH * diff;
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


} // namespace cg::snd