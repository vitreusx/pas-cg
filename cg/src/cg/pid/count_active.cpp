#include <cg/pid/count_active.h>

namespace cg::pid {
active_counts count_active::operator()() const {
  active_counts counts;

  for (auto const &bundle : *eval->bundles) {
    auto i1 = bundle.i1(), i2 = bundle.i2();
    auto chain1 = chain_idx[i1], chain2 = chain_idx[i2];
    vec3r r1 = eval->r[i1], r2 = eval->r[i2];
    auto r12 = eval->simul_box->wrap(r1, r2);
    auto r12_rn = norm_inv(r12);
    if ((real)1.0 > r12_rn * eval->cutoff)
      continue;

    real psi1, psi2;
    vec3r dpsi1_dr1p, dpsi1_dr1, dpsi1_dr1n, dpsi1_dr2;
    vec3r dpsi2_dr2p, dpsi2_dr2, dpsi2_dr2n, dpsi2_dr1;

    auto type = bundle.type();
    //  auto i1p = prev[i1], i1n = next[i1], i2p = prev[i2], i2n = next[i2];
    //  if (i1p < 0 || i1n < 0 || i2p < 0 || i2n < 0)
    //    return;
    auto i1p = i1 - 1, i1n = i1 + 1, i2p = i2 - 1, i2n = i2 + 1;
    auto n = eval->r.size();

    lambda const *bb_lam_1, *bb_lam_2;

    {
      auto i1_ = i1, i2_ = i1n, i3_ = i1p, i4_ = i2;
      //    auto r1_ = r[i1_], r2_ = r[i2_], r3_ = r[i3_], r4_ = r[i4_];
      //    auto rij = r1_ - r2_, rkj = r3_ - r2_, rkl = r3_ - r4_;
      auto rij = eval->r[i1_] - (i2_ < n ? eval->r[i2_] : vec3r());
      auto rkj = i3_ >= 0 ? eval->r[i3_] - (i2_ < n ? eval->r[i2_] : vec3r())
                          : vec3r();
      auto rkl = (i3_ >= 0 ? eval->r[i3_] : vec3r()) - eval->r[i4_];
      auto rm = cross(rij, rkj);
      auto rn = cross(rkj, rkl);
      auto rm_ninv = norm_inv(rm), rn_ninv = norm_inv(rn),
           rkj_ninv = norm_inv(rkj);
      if (rm_ninv > (real)10.0 || rn_ninv > (real)10.0)
        continue;
      auto rkj_n = norm(rkj);
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
        bb_lam_1 = &eval->bb_plus_lam;
      } else {
        bb_lam_1 = &eval->bb_minus_lam;
      }
      if (!bb_lam_1->supp(psi1) && (i2 - i1 != 3 || m != 1)) {
        bb_lam_1 = &eval->bb_minus_lam;
      }
    }

    {
      auto i1_ = i2, i2_ = i2n, i3_ = i2p, i4_ = i1;
      //    auto r1_ = r[i1_], r2_ = r[i2_], r3_ = r[i3_], r4_ = r[i4_];
      //    auto rij = r1_ - r2_, rkj = r3_ - r2_, rkl = r3_ - r4_;
      auto rij = eval->r[i1_] - (i2_ < n ? eval->r[i2_] : vec3r());
      auto rkj = i3_ >= 0 ? eval->r[i3_] - (i2_ < n ? eval->r[i2_] : vec3r())
                          : vec3r();
      auto rkl = (i3_ >= 0 ? eval->r[i3_] : vec3r()) - eval->r[i4_];
      auto rm = cross(rij, rkj);
      auto rn = cross(rkj, rkl);
      auto rm_ninv = norm_inv(rm), rn_ninv = norm_inv(rn),
           rkj_ninv = norm_inv(rkj);
      if (rm_ninv > (real)10.0 || rn_ninv > (real)10.0)
        continue;
      auto rkj_n = norm(rkj);
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
        bb_lam_2 = &eval->bb_plus_lam;
      } else {
        bb_lam_2 = &eval->bb_minus_lam;
      }
      if (!bb_lam_2->supp(psi2) && (i2 - i1 != 3 || m != 1)) {
        bb_lam_2 = &eval->bb_minus_lam;
      }
    }

    if (bb_lam_1->supp(psi1) && bb_lam_2->supp(psi2)) {
      auto [lam1, deriv1] = (*bb_lam_1)(psi1);
      auto [lam2, deriv2] = (*bb_lam_2)(psi2);

      if (lam1 * lam2 > (real)5e-5) {
        if (chain1 == chain2) {
          ++counts.intra;
          ++counts.intra_bb;
        }
        else {
          ++counts.inter;
          ++counts.inter_bb;
        }
      }
    }

    sink_lj ss_sink_lj = eval->ss_ljs[type];
    if (eval->ss_lam.supp(psi1) && eval->ss_lam.supp(psi2)) {
      auto [lam1, deriv1] = eval->ss_lam(psi1);
      auto [lam2, deriv2] = eval->ss_lam(psi2);

      if (lam1 * lam2 > (real)5e-5 && ss_sink_lj.depth() > (real)5e-5) {
        if (chain1 == chain2) {
          ++counts.intra;
          ++counts.intra_ss;
        }
        else {
          ++counts.inter;
          ++counts.inter_ss;
        }
      }
    }
  }

  return counts;
}
} // namespace cg::pid