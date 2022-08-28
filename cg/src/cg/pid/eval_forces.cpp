#include <cg/pid/eval_forces.h>
#include <iostream>

namespace cg::pid {

void eval_forces::operator()() {
  for (int idx = 0; idx < bundles->size(); ++idx) {
    iter(bundles->at(idx));
  }
}

template <typename E>
void eval_forces::iter(bundle_expr<E> const &bundle) const {
  if (bundle.orig_dist() - *total_disp > cutoff)
    return;

  auto i1 = bundle.i1(), i2 = bundle.i2();
  if ((aa_code)atype[i1] == aa_code::PRO || (aa_code)atype[i2] == aa_code::PRO)
    return;

  vec3r r1 = r[i1], r2 = r[i2];
  auto r12 = simul_box->wrap(r1, r2);
  auto r12_rn = norm_inv(r12);
  if ((real)1.0 > r12_rn * cutoff)
    return;

  auto r12_n = (real)1.0 / r12_rn;

  real psi1, psi2;
  vec3r dpsi1_dr1p, dpsi1_dr1, dpsi1_dr1n, dpsi1_dr2;
  vec3r dpsi2_dr2p, dpsi2_dr2, dpsi2_dr2n, dpsi2_dr1;

  auto type = bundle.type();
  //  auto i1p = prev[i1], i1n = next[i1], i2p = prev[i2], i2n = next[i2];
  //  if (i1p < 0 || i1n < 0 || i2p < 0 || i2n < 0)
  //    return;
  auto i1p = i1 - 1, i1n = i1 + 1, i2p = i2 - 1, i2n = i2 + 1;
  auto n = r.size();

  lambda const *bb_lam_1, *bb_lam_2;
  sink_lj const *bb_lj;

  {
    auto i1_ = i1, i2_ = i1n, i3_ = i1p, i4_ = i2;
    //    auto r1_ = r[i1_], r2_ = r[i2_], r3_ = r[i3_], r4_ = r[i4_];
    //    auto rij = r1_ - r2_, rkj = r3_ - r2_, rkl = r3_ - r4_;
    auto rij = r[i1_] - (i2_ < n ? r[i2_] : vec3r());
    auto rkj = i3_ >= 0 ? r[i3_] - (i2_ < n ? r[i2_] : vec3r()) : vec3r();
    auto rkl = (i3_ >= 0 ? r[i3_] : vec3r()) - r[i4_];

    auto rm = cross(rij, rkj);
    auto rn = cross(rkj, rkl);

    auto rm_n = norm(rm), rn_n = norm(rn);
    if (rm_n < (real)0.1 || rn_n < (real)0.1)
      return;

    auto rm_ninv = (real)1.0 / rm_n, rn_ninv = (real)1.0 / rn_n,
         rkj_ninv = norm_inv(rkj), rkj_n = (real)1.0 / rkj_ninv;

    auto fi = rm * rkj_n * rm_ninv * rm_ninv;
    auto fl = -rn * rkj_n * rn_ninv * rn_ninv;
    auto df = (fi * dot(rij, rkj) - fl * dot(rkl, rkj)) * rkj_ninv * rkj_ninv;
    auto fj = -fi + df;
    auto fk = -fl - df;

    dpsi1_dr1 = fi;
    dpsi1_dr1n = fj;
    dpsi1_dr1p = fk;
    dpsi1_dr2 = fl;

    auto cos_psi1 = dot(rm, rn) * rm_ninv * rn_ninv;
    psi1 = acos(cos_psi1);
    if (dot(rij, rn) < (real)0.0)
      psi1 = -psi1;

    int m = 1;
    if (i2 - i1 != 3 || m != 2) {
      bb_lam_1 = &bb_plus_lam;
      bb_lj = &bb_plus_lj;
    } else {
      bb_lam_1 = &bb_minus_lam;
      bb_lj = &bb_minus_lj;
    }
    if (!bb_lam_1->supp(psi1) && (i2 - i1 != 3 || m != 1)) {
      bb_lam_1 = &bb_minus_lam;
      bb_lj = &bb_minus_lj;
    }
  }

  {
    auto i1_ = i2, i2_ = i2n, i3_ = i2p, i4_ = i1;
    //    auto r1_ = r[i1_], r2_ = r[i2_], r3_ = r[i3_], r4_ = r[i4_];
    //    auto rij = r1_ - r2_, rkj = r3_ - r2_, rkl = r3_ - r4_;
    auto rij = r[i1_] - (i2_ < n ? r[i2_] : vec3r());
    auto rkj = i3_ >= 0 ? r[i3_] - (i2_ < n ? r[i2_] : vec3r()) : vec3r();
    auto rkl = (i3_ >= 0 ? r[i3_] : vec3r()) - r[i4_];

    auto rm = cross(rij, rkj);
    auto rn = cross(rkj, rkl);

    auto rm_n = norm(rm), rn_n = norm(rn);
    if (rm_n < (real)0.1 || rn_n < (real)0.1)
      return;

    auto rm_ninv = (real)1.0 / rm_n, rn_ninv = (real)1.0 / rn_n,
         rkj_ninv = norm_inv(rkj), rkj_n = (real)1.0 / rkj_ninv;

    auto fi = rm * rkj_n * rm_ninv * rm_ninv;
    auto fl = -rn * rkj_n * rn_ninv * rn_ninv;
    auto df = (fi * dot(rij, rkj) - fl * dot(rkl, rkj)) * rkj_ninv * rkj_ninv;
    auto fj = -fi + df;
    auto fk = -fl - df;

    dpsi2_dr2 = fi;
    dpsi2_dr2n = fj;
    dpsi2_dr2p = fk;
    dpsi2_dr1 = fl;

    auto cos_psi2 = dot(rm, rn) * rm_ninv * rn_ninv;
    psi2 = acos(cos_psi2);
    if (dot(rij, rn) < (real)0.0)
      psi2 = -psi2;

    int m = 2;
    if (i2 - i1 != 3 || m != 2) {
      bb_lam_2 = &bb_plus_lam;
      bb_lj = &bb_plus_lj;
    } else {
      bb_lam_2 = &bb_minus_lam;
      bb_lj = &bb_minus_lj;
    }
    if (!bb_lam_2->supp(psi2) && (i2 - i1 != 3 || m != 1)) {
      bb_lam_2 = &bb_minus_lam;
      bb_lj = &bb_minus_lj;
    }
  }

  real dV_dpsi1 = 0.0f, dV_dpsi2 = 0.0f, dV_dr = 0.0f, V_ = 0.0f;

  if (bb_lam_1->supp(psi1) && bb_lam_2->supp(psi2)) {
    auto [lam1, deriv1] = (*bb_lam_1)(psi1);
    auto [lam2, deriv2] = (*bb_lam_2)(psi2);

    if (lam1 * lam2 > (real)5e-5) {
      auto [lj_V, lj_dV_dr] = (*bb_lj)(r12_n, r12_rn);

      V_ += lam1 * lam2 * lj_V;
      dV_dpsi1 += deriv1 * lam2 * lj_V;
      dV_dpsi2 += deriv2 * lam1 * lj_V;
      dV_dr += lam1 * lam2 * lj_dV_dr;
    }
  }

  sink_lj ss_sink_lj = ss_ljs[type];
  if (ss_lam.supp(psi1) && ss_lam.supp(psi2)) {
    auto [lam1, deriv1] = ss_lam(psi1);
    auto [lam2, deriv2] = ss_lam(psi2);

    if (lam1 * lam2 > (real)5e-5 && ss_sink_lj.depth() > (real)5e-5) {
      auto [lj_V, lj_dV_dr] = ss_sink_lj(r12_n, r12_rn);

      V_ += lam1 * lam2 * lj_V;
      dV_dpsi1 += deriv1 * lam2 * lj_V;
      dV_dpsi2 += deriv2 * lam1 * lj_V;
      dV_dr += lam1 * lam2 * lj_dV_dr;
    }
  }

  *V += V_;
  //  if (dV_dr != 0.0) {
  //    auto fmt = std::cout.flags();
  //    std::cout << std::scientific << i1 + 1 << " " << i2 + 1 << " " << dV_dr
  //              << '\n';
  //    std::cout.flags(fmt);
  //  }

  auto r12_u = r12 * r12_rn;
  if (i1p >= 0)
    F[i1p] -= dV_dpsi1 * dpsi1_dr1p;
  F[i1] -= dV_dpsi1 * dpsi1_dr1 + dV_dpsi2 * dpsi2_dr1 - dV_dr * r12_u;
  if (i1n < n)
    F[i1n] -= dV_dpsi1 * dpsi1_dr1n;
  if (i2p >= 0)
    F[i2p] -= dV_dpsi2 * dpsi2_dr2p;
  F[i2] -= dV_dpsi1 * dpsi1_dr2 + dV_dpsi2 * dpsi2_dr2 + dV_dr * r12_u;
  if (i2n < n)
    F[i2n] -= dV_dpsi2 * dpsi2_dr2n;
}

void eval_forces::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < bundles->size(); ++idx) {
    iter(bundles->at(idx));
  }
}

bool eval_forces::is_active(const bundle &bundle) const {
  auto i1 = bundle.i1(), i2 = bundle.i2();
  vec3r r1 = r[i1], r2 = r[i2];

  auto r12 = norm(simul_box->wrap(r1, r2));
  auto r_min = ss_ljs[bundle.type()].r_high();
  auto sigma = r_min * C216_INV;
  return r12 <= active_thr * sigma;
}

void eval_forces::for_slice(int from, int to) const {
  for (int idx = from; idx < to; ++idx)
    iter(bundles->at(idx));
}

int eval_forces::total_size() const {
  return bundles->size();
}


} // namespace cg::pid