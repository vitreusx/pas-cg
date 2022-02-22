#include "features/pid/eval_forces.h"
using namespace cg::pid;

void eval_forces::operator()() {
  for (int idx = 0; idx < bundles->size(); ++idx) {
    iter(bundles->at(idx));
  }
}

template <typename E>
void eval_forces::iter(bundle_expr<E> const &bundle) const {
  real psi1, psi2;
  vec3r dpsi1_dr1p, dpsi1_dr1, dpsi1_dr1n, dpsi1_dr2;
  vec3r dpsi2_dr2p, dpsi2_dr2, dpsi2_dr2n, dpsi2_dr1;

  auto i1 = bundle.i1(), i2 = bundle.i2();
  auto type = bundle.type();
  auto i1p = prev->at(i1), i1n = next->at(i1), i2p = prev->at(i2),
       i2n = next->at(i2);

  vec3r r1p = r->at(i1p), r1 = r->at(i1), r1n = r->at(i1n);
  vec3r r2p = r->at(i2p), r2 = r->at(i2), r2n = r->at(i2n);

  {
    auto &r1_ = r1p, &r2_ = r1, &r3_ = r1n, &r4_ = r2;
    auto r24 = box->r_uv(r2_, r4_);
    auto r12 = r2_ - r1_;
    auto r23 = r3_ - r2_;
    auto r13 = r3_ - r1_;
    auto r14 = r12 + r24;

    auto rij = -r23;
    auto rkj = -r13;
    auto rkl = -r14;
    auto rm = cross(r23, r12);
    auto rn = cross(r13, r14);
    auto rm_rn = norm_inv(rm);
    auto rn_rn = norm_inv(rn);
    auto rkj_n = norm(rkj);

    auto fi = rm * rkj_n * rm_rn;
    auto fl = -rn * rkj_n * rn_rn;
    auto df = (fi * dot(rij, rkj) - fl * dot(rkl, rkj)) / (rkj_n * rkj_n);
    auto fj = -fi + df;
    auto fk = -fl - df;

    dpsi1_dr1p = fi;
    dpsi1_dr1 = fj;
    dpsi1_dr1n = fk;
    dpsi1_dr2 = fl;

    auto cos_psi1 = dot(rm, rn) * rm_rn * rn_rn;
    psi1 = acos(cos_psi1);
    if (dot(rij, rn) < 0.0f)
      psi1 = -psi1;
  }

  {
    auto &r1_ = r2p, &r2_ = r2, &r3_ = r2n, &r4_ = r1;
    auto r24 = box->r_uv(r2_, r4_);
    auto r12 = r2_ - r1_;
    auto r23 = r3_ - r2_;
    auto r13 = r3_ - r1_;
    auto r14 = r12 + r24;

    auto rij = -r23;
    auto rkj = -r13;
    auto rkl = -r14;
    auto rm = cross(r23, r12);
    auto rn = cross(r13, r14);
    auto rm_rn = norm_inv(rm);
    auto rn_rn = norm_inv(rn);
    auto rkj_n = norm(rkj);

    auto fi = rm * rkj_n * rm_rn;
    auto fl = -rn * rkj_n * rn_rn;
    auto df = (fi * dot(rij, rkj) - fl * dot(rkl, rkj)) / (rkj_n * rkj_n);
    auto fj = -fi + df;
    auto fk = -fl - df;

    dpsi1_dr1p = fi;
    dpsi1_dr1 = fj;
    dpsi1_dr1n = fk;
    dpsi1_dr2 = fl;

    auto cos_psi2 = dot(rm, rn) * rm_rn * rn_rn;
    psi2 = acos(cos_psi2);
    if (dot(rij, rn) < 0.0f)
      psi2 = -psi2;
  }

  auto r12 = box->r_uv(r1, r2);
  auto r12_n = norm(r12), r12_rn = 1.0f / r12_n;

  real A = 0.0f, B = 0.0f, C = 0.0f, V_ = 0.0f;

  if (bb_plus_lam.supp(psi1) && bb_plus_lam.supp(psi2)) {
    auto [lam1, deriv1] = bb_plus_lam(psi1);
    auto [lam2, deriv2] = bb_plus_lam(psi2);
    auto [lj_V, lj_dV_dr] = bb_plus_lj(r12_rn);

    V_ += lam1 * lam2 * lj_V;
    A += deriv1 * lam2 * lj_V;
    B += deriv2 * lam1 * lj_V;
    C += lam1 * lam2 * lj_dV_dr;
  }

  if (bb_minus_lam.supp(psi1) && bb_minus_lam.supp(psi2)) {
    auto [lam1, deriv1] = bb_minus_lam(psi1);
    auto [lam2, deriv2] = bb_minus_lam(psi2);
    auto [lj_V, lj_dV_dr] = bb_minus_lj(r12_rn);

    V_ += lam1 * lam2 * lj_V;
    A += deriv1 * lam2 * lj_V;
    B += deriv2 * lam1 * lj_V;
    C += lam1 * lam2 * lj_dV_dr;
  }

  sink_lj ss_sink_lj = ss_ljs->at(type);
  if (ss_lam.supp(psi1) && ss_lam.supp(psi2)) {
    auto [lam1, deriv1] = ss_lam(psi1);
    auto [lam2, deriv2] = ss_lam(psi2);
    auto [lj_V, lj_dV_dr] = ss_sink_lj(r12_n, r12_rn);

    V_ += lam1 * lam2 * lj_V;
    A += deriv1 * lam2 * lj_V;
    B += deriv2 * lam1 * lj_V;
    C += lam1 * lam2 * lj_dV_dr;
  }

  *V += V_;

  auto r12_u = r12 * r12_rn;
  F->at(i1p) -= A * dpsi1_dr1p;
  F->at(i1) -= A * dpsi1_dr1 + B * dpsi2_dr1 + C * r12_u;
  F->at(i1n) -= A * dpsi1_dr1n;
  F->at(i2p) -= B * dpsi2_dr2p;
  F->at(i2) -= A * dpsi1_dr2 + B * dpsi2_dr2 - C * r12_u;
  F->at(i2n) -= B * dpsi2_dr2n;
}

void eval_forces::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < bundles->size(); ++idx) {
    iter(bundles->at(idx));
  }
}