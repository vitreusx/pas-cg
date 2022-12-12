#include <cg/pid/eval_forces.h>
#include <cg/pid/vect_utils.h>

namespace cg::pid {
void eval_forces::iter(int idx) const {
  //  if (bundle.orig_dist() - *total_disp > cutoff)
  //    return;

  auto bundle = (*bundles)[idx];
  auto i1 = bundle.i1(), i2 = bundle.i2();
  if ((aa_code)atype[i1] == aa_code::PRO || (aa_code)atype[i2] == aa_code::PRO)
    return;

  vec3r r1 = r[i1], r2 = r[i2];
  auto r12 = simul_box->wrap<vec3r>(r1, r2);
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
    auto rij = r[i1_] - (i2_ < n ? (vec3r)r[i2_] : vec3r());
    auto rkj =
        i3_ >= 0 ? r[i3_] - (i2_ < n ? (vec3r)r[i2_] : vec3r()) : vec3r();
    auto rkl = (i3_ >= 0 ? (vec3r)r[i3_] : vec3r()) - r[i4_];

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
    auto rij = r[i1_] - (i2_ < n ? (vec3r)r[i2_] : vec3r());
    auto rkj =
        i3_ >= 0 ? r[i3_] - (i2_ < n ? (vec3r)r[i2_] : vec3r()) : vec3r();
    auto rkl = (i3_ >= 0 ? (vec3r)r[i3_] : vec3r()) - r[i4_];

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

void eval_forces::vect_iter(int lane_idx) const {
  static constexpr int N = elems_per_vect, W = vect::VECT_BITS;

  using mask_t = vect::lane<bool, N, W>;
  using reals = vect::lane<real, N, W>;
  using vec3rs = vect::lane<vec3r, N, W>;

  auto bundle = bundles->template at_lane<N, W>(lane_idx);
  auto i1 = bundle.i1(), i2 = bundle.i2();

  auto mask =
      mask_t(atype[i1] != aa_code::PRO) & mask_t(atype[i2] != aa_code::PRO);

  auto r1 = r[i1], r2 = r[i2];
  auto r12 = simul_box->wrap<vect::lane<vec3r, N, W>>(r1, r2);
  auto r12_n = norm_(r12);
  auto r12_rn = (real)1.0 / r12_n;
  mask = mask & mask_t(r12_n < cutoff);

  vect::lane<real, N, W> psi1, psi2;
  vect::lane<vec3r, N, W> dpsi1_dr1p, dpsi1_dr1, dpsi1_dr1n, dpsi1_dr2;
  vect::lane<vec3r, N, W> dpsi2_dr2p, dpsi2_dr2, dpsi2_dr2n, dpsi2_dr1;

  vect::lane<int, N, W> i1p = i1 - 1, i1n = i1 + 1, i2p = i2 - 1, i2n = i2 + 1;

  auto n = r.size();
  auto val_i1p = (i1p >= 0), val_i1n = (i1n < n), val_i2p = (i2p >= 0),
       val_i2n = (i2n < n);

  auto r_i1p = r[std::make_pair(i1p, val_i1p)],
       r_i1n = r[std::make_pair(i1n, val_i1n)],
       r_i2p = r[std::make_pair(i2p, val_i2p)],
       r_i2n = r[std::make_pair(i2n, val_i2n)];

  vect::lane<bool, N, W> bb_lam_1_opt, bb_lam_2_opt, bb_lj_opt;

  vect::lane<lambda, N, W> bb_plus_lam_v = bb_plus_lam,
                           bb_minus_lam_v = bb_minus_lam;

  {
    //    auto &i1_ = i1, &i2_ = i1n, &i3_ = i1p, &i4_ = i2;
    //    auto rkj =
    //        select(i3_ >= 0, r[i3_] - select(i2_ < n, r[i2_], zero_v),
    //        zero_v);
    //    auto rkl = select(i3_ >= 0, r[i3_], zero_v) - r[i4_];

    auto rij = r1 - r_i1n;
    auto rkj = select_(val_i1p, vec3rs(r_i1p - r_i1n), vec3rs());
    auto rkl = r_i1p - r2;

    auto rm = cross(rij, rkj);
    auto rn = cross(rkj, rkl);

    auto rm_n = norm_(rm), rn_n = norm_(rn);
    mask = mask & mask_t(rm_n >= (real)0.1) & mask_t(rn_n >= (real)0.1);

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
    cos_psi1 = ::max(::min(cos_psi1, (reals)1.0), (reals)-1.0);
    psi1 = acos(cos_psi1);
    psi1 = select_(dot(rij, rn) < (reals)0.0, -psi1, psi1);

    int m = 1;
    bb_lam_1_opt = (i2 - i1 != 3) | mask_t(m != 2);

    auto bb_lam_1 = select_(bb_lam_1_opt, bb_plus_lam_v, bb_minus_lam_v);
    bb_lam_1_opt = bb_lam_1_opt & !((!mask_t(bb_lam_1.supp(psi1))) &
                                    mask_t((i2 - i1 != 3) | mask_t(m != 1)));
  }

  {
    //    auto &i1_ = i2, &i2_ = i2n, &i3_ = i2p, &i4_ = i1;
    //    auto rij = r[i1_] - select(i2_ < n, r[i2_], zero_v);
    //    auto rkj =
    //        select(i3_ >= 0, r[i3_] - select(i2_ < n, r[i2_], zero_v),
    //        zero_v);
    //    auto rkl = select(i3_ >= 0, r[i3_], zero_v) - r[i4_];

    auto rij = r2 - r_i2n;
    auto rkj = select_(val_i2p, vec3rs(r_i2p - r_i2n), vec3rs());
    auto rkl = r_i2p - r1;

    auto rm = cross(rij, rkj);
    auto rn = cross(rkj, rkl);

    auto rm_n = norm_(rm), rn_n = norm_(rn);
    mask = mask & mask_t(rm_n >= (real)0.1) & mask_t(rn_n >= (real)0.1);

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
    cos_psi2 = ::max(::min(cos_psi2, (reals)1.0), (reals)-1.0);
    psi2 = acos(cos_psi2);
    psi2 = select_(dot(rij, rn) < (reals)0.0, -psi2, psi2);

    int m = 2;
    bb_lam_2_opt = (i2 - i1 != 3) | mask_t(m != 2);

    auto bb_lam_2 = select_(bb_lam_2_opt, bb_plus_lam_v, bb_minus_lam_v);
    bb_lam_2_opt = bb_lam_2_opt & !((!mask_t(bb_lam_2.supp(psi2))) &
                                    mask_t((i2 - i1 != 3) | mask_t(m != 1)));

    bb_lj_opt = bb_lam_2_opt;
  }

  vect::lane<real, N, W> dV_dpsi1 = (real)0.0, dV_dpsi2 = (real)0.0,
                         dV_dr = (real)0.0, V_ = (real)0.0;

  auto bb_lam_1 = select_(bb_lam_1_opt, bb_plus_lam_v, bb_minus_lam_v);
  auto bb_lam_2 = select_(bb_lam_2_opt, bb_plus_lam_v, bb_minus_lam_v);

  mask_t bb_supp = (mask_t)(bb_lam_1.supp(psi1) & bb_lam_2.supp(psi2));
  if (horizontal_or(bb_supp)) {
    auto [lam1, deriv1] = bb_lam_1(psi1);
    auto [lam2, deriv2] = bb_lam_2(psi2);

    vect::lane<sink_lj, N, W> bb_minus_lj_v = bb_minus_lj,
                              bb_plus_lj_v = bb_plus_lj;
    auto bb_lj = select_(bb_lj_opt, bb_plus_lj_v, bb_minus_lj_v);

    auto _mask = bb_supp & mask_t(lam1 * lam2 > (real)5e-5);
    auto [lj_V, lj_dV_dr] = bb_lj(r12_n, r12_rn);

    V_ = if_add(_mask, V_, lam1 * lam2 * lj_V);
    dV_dpsi1 = if_add(_mask, dV_dpsi1, deriv1 * lam2 * lj_V);
    dV_dpsi2 = if_add(_mask, dV_dpsi2, deriv2 * lam1 * lj_V);
    dV_dr = if_add(_mask, dV_dr, lam1 * lam2 * lj_dV_dr);
  }

  auto type = bundle.type();
  auto ss_sink_lj = ss_ljs[type];

  mask_t ss_supp = (mask_t)(ss_lam.supp(psi1) & ss_lam.supp(psi2));
  if (horizontal_or(ss_supp)) {
    auto [lam1, deriv1] = ss_lam(psi1);
    auto [lam2, deriv2] = ss_lam(psi2);

    auto _mask = ss_supp & mask_t(lam1 * lam2 > (real)5e-5) &
                 mask_t(ss_sink_lj.depth() > (real)5e-5);
    auto [lj_V, lj_dV_dr] = ss_sink_lj(r12_n, r12_rn);

    V_ = if_add(_mask, V_, lam1 * lam2 * lj_V);
    dV_dpsi1 = if_add(_mask, dV_dpsi1, deriv1 * lam2 * lj_V);
    dV_dpsi2 = if_add(_mask, dV_dpsi2, deriv2 * lam1 * lj_V);
    dV_dr = if_add(_mask, dV_dr, lam1 * lam2 * lj_dV_dr);
  }

  auto r12_u = r12 * r12_rn;

  vect::lane<real, N, W> zero_v = (real)0;
  V_ = select_(mask, V_, zero_v);
  *V += horizontal_add(V_);

  auto diff_F_i1p = to_array<vec3r, N, W>(dV_dpsi1 * dpsi1_dr1p);
  auto diff_F_i1 = to_array<vec3r, N, W>(dV_dpsi1 * dpsi1_dr1 +
                                         dV_dpsi2 * dpsi2_dr1 - dV_dr * r12_u);
  auto diff_F_i1n = to_array<vec3r, N, W>(dV_dpsi1 * dpsi1_dr1n);
  auto diff_F_i2p = to_array<vec3r, N, W>(dV_dpsi2 * dpsi2_dr2p);
  auto diff_F_i2 = to_array<vec3r, N, W>(dV_dpsi1 * dpsi1_dr2 +
                                         dV_dpsi2 * dpsi2_dr2 + dV_dr * r12_u);
  auto diff_F_i2n = to_array<vec3r, N, W>(dV_dpsi2 * dpsi2_dr2n);

  for (int idx = 0; idx < (int)N; ++idx) {
    if (!mask[idx])
      continue;

    if (val_i1p[idx])
      F[i1p[idx]] -= diff_F_i1p[idx];
    F[i1[idx]] -= diff_F_i1[idx];
    if (val_i1n[idx])
      F[i1n[idx]] -= diff_F_i1n[idx];
    if (val_i2p[idx])
      F[i2p[idx]] -= diff_F_i2p[idx];
    F[i2[idx]] -= diff_F_i2[idx];
    if (val_i2n[idx])
      F[i2n[idx]] -= diff_F_i2n[idx];
  }
}

int eval_forces::size() const {
  return bundles->size();
}

bool eval_forces::is_active(const bundle &bundle) const {
  auto i1 = bundle.i1(), i2 = bundle.i2();
  vec3r r1 = r[i1], r2 = r[i2];

  auto r12 = norm(simul_box->wrap<vec3r>(r1, r2));
  auto r_min = ss_ljs[bundle.type()].r_high();
  auto sigma = r_min * C216_INV;
  return r12 <= active_thr * sigma;
}

void eval_forces::fast_vect_iter(int lane_idx) const {
  static constexpr int N = elems_per_vect, W = vect::VECT_BITS;

  using mask_t = vect::lane<bool, N, W>;
  using reals = vect::lane<real, N, W>;

  auto bundle = bundles->template at_lane<N, W>(lane_idx);
  auto i1 = bundle.i1(), i2 = bundle.i2();

  auto mask =
      mask_t(atype[i1] != aa_code::PRO) & mask_t(atype[i2] != aa_code::PRO);

  auto r1 = r[i1], r2 = r[i2];
  auto r12 = simul_box->wrap<vect::lane<vec3r, N, W>>(r1, r2);
  auto r12_n = norm_(r12);
  auto r12_rn = (real)1.0 / r12_n;
  mask = mask & mask_t(r12_n < cutoff);

  vect::lane<real, N, W> psi1, psi2;
  vect::lane<vec3r, N, W> dpsi1_dr1p, dpsi1_dr1, dpsi1_dr1n, dpsi1_dr2;
  vect::lane<vec3r, N, W> dpsi2_dr2p, dpsi2_dr2, dpsi2_dr2n, dpsi2_dr1;

  vect::lane<int, N, W> i1p = i1 - 1, i1n = i1 + 1, i2p = i2 - 1, i2n = i2 + 1;

  auto r_i1p = r[i1p], r_i1n = r[i1n], r_i2p = r[i2p], r_i2n = r[i2n];

  vect::lane<bool, N, W> bb_lam_1_opt, bb_lam_2_opt, bb_lj_opt;

  vect::lane<lambda, N, W> bb_plus_lam_v = bb_plus_lam,
                           bb_minus_lam_v = bb_minus_lam;

  {
    //    auto &i1_ = i1, &i2_ = i1n, &i3_ = i1p, &i4_ = i2;
    //    auto rkj =
    //        select(i3_ >= 0, r[i3_] - select(i2_ < n, r[i2_], zero_v),
    //        zero_v);
    //    auto rkl = select(i3_ >= 0, r[i3_], zero_v) - r[i4_];

    auto rij = r1 - r_i1n;
    auto rkj = r_i1p - r_i1n;
    auto rkl = r_i1p - r2;

    auto rm = cross(rij, rkj);
    auto rn = cross(rkj, rkl);

    auto rm_n = norm_(rm), rn_n = norm_(rn);
    mask = mask & mask_t(rm_n >= (real)0.1) & mask_t(rn_n >= (real)0.1);

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
    cos_psi1 = ::max(::min(cos_psi1, (reals)1.0), (reals)-1.0);
    psi1 = acos(cos_psi1);
    psi1 = select_(dot(rij, rn) < (reals)0.0, -psi1, psi1);

    int m = 1;
    bb_lam_1_opt = (i2 - i1 != 3) | mask_t(m != 2);

    auto bb_lam_1 = select_(bb_lam_1_opt, bb_plus_lam_v, bb_minus_lam_v);
    bb_lam_1_opt = bb_lam_1_opt & !((!mask_t(bb_lam_1.supp(psi1))) &
                                    mask_t((i2 - i1 != 3) | mask_t(m != 1)));
  }

  {
    //    auto &i1_ = i2, &i2_ = i2n, &i3_ = i2p, &i4_ = i1;
    //    auto rij = r[i1_] - select(i2_ < n, r[i2_], zero_v);
    //    auto rkj =
    //        select(i3_ >= 0, r[i3_] - select(i2_ < n, r[i2_], zero_v),
    //        zero_v);
    //    auto rkl = select(i3_ >= 0, r[i3_], zero_v) - r[i4_];

    auto rij = r2 - r_i2n;
    auto rkj = r_i2p - r_i2n;
    auto rkl = r_i2p - r1;

    auto rm = cross(rij, rkj);
    auto rn = cross(rkj, rkl);

    auto rm_n = norm_(rm), rn_n = norm_(rn);
    mask = mask & mask_t(rm_n >= (real)0.1) & mask_t(rn_n >= (real)0.1);

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
    cos_psi2 = ::max(::min(cos_psi2, (reals)1.0), (reals)-1.0);
    psi2 = acos(cos_psi2);
    psi2 = select_(dot(rij, rn) < (reals)0.0, -psi2, psi2);

    int m = 2;
    bb_lam_2_opt = (i2 - i1 != 3) | mask_t(m != 2);

    auto bb_lam_2 = select_(bb_lam_2_opt, bb_plus_lam_v, bb_minus_lam_v);
    bb_lam_2_opt = bb_lam_2_opt & !((!mask_t(bb_lam_2.supp(psi2))) &
                                    mask_t((i2 - i1 != 3) | mask_t(m != 1)));

    bb_lj_opt = bb_lam_2_opt;
  }

  vect::lane<real, N, W> dV_dpsi1 = (real)0.0, dV_dpsi2 = (real)0.0,
                         dV_dr = (real)0.0, V_ = (real)0.0;

  auto bb_lam_1 = select_(bb_lam_1_opt, bb_plus_lam_v, bb_minus_lam_v);
  auto bb_lam_2 = select_(bb_lam_2_opt, bb_plus_lam_v, bb_minus_lam_v);

  mask_t bb_supp = (mask_t)(mask & bb_lam_1.supp(psi1) & bb_lam_2.supp(psi2));
  if (horizontal_or(bb_supp)) {
    auto [lam1, deriv1] = bb_lam_1(psi1);
    auto [lam2, deriv2] = bb_lam_2(psi2);

    vect::lane<sink_lj, N, W> bb_minus_lj_v = bb_minus_lj,
                              bb_plus_lj_v = bb_plus_lj;
    auto bb_lj = select_(bb_lj_opt, bb_plus_lj_v, bb_minus_lj_v);

    auto _mask = bb_supp & mask_t(lam1 * lam2 > (real)5e-5);
    auto [lj_V, lj_dV_dr] = bb_lj(r12_n, r12_rn);

    V_ = if_add(_mask, V_, lam1 * lam2 * lj_V);
    dV_dpsi1 = if_add(_mask, dV_dpsi1, deriv1 * lam2 * lj_V);
    dV_dpsi2 = if_add(_mask, dV_dpsi2, deriv2 * lam1 * lj_V);
    dV_dr = if_add(_mask, dV_dr, lam1 * lam2 * lj_dV_dr);
  }

  auto type = bundle.type();
  auto ss_sink_lj = ss_ljs[type];

  mask_t ss_supp = (mask_t)(mask & ss_lam.supp(psi1) & ss_lam.supp(psi2));
  if (horizontal_or(ss_supp)) {
    auto [lam1, deriv1] = ss_lam(psi1);
    auto [lam2, deriv2] = ss_lam(psi2);

    auto _mask = ss_supp & mask_t(lam1 * lam2 > (real)5e-5) &
                 mask_t(ss_sink_lj.depth() > (real)5e-5);
    auto [lj_V, lj_dV_dr] = ss_sink_lj(r12_n, r12_rn);

    V_ = if_add(_mask, V_, lam1 * lam2 * lj_V);
    dV_dpsi1 = if_add(_mask, dV_dpsi1, deriv1 * lam2 * lj_V);
    dV_dpsi2 = if_add(_mask, dV_dpsi2, deriv2 * lam1 * lj_V);
    dV_dr = if_add(_mask, dV_dr, lam1 * lam2 * lj_dV_dr);
  }

  auto r12_u = r12 * r12_rn;

  vect::lane<real, N, W> zero_v = (real)0;
  V_ = select_(mask, V_, zero_v);
  *V += horizontal_add(V_);

  F[i1p] = F[i1p] - dV_dpsi1 * dpsi1_dr1p;
  F[i1] = F[i1] - (dV_dpsi1 * dpsi1_dr1 + dV_dpsi2 * dpsi2_dr1 - dV_dr * r12_u);
  F[i1n] = F[i1n] - dV_dpsi1 * dpsi1_dr1n;
  F[i2p] = F[i2p] - dV_dpsi2 * dpsi2_dr2p;
  F[i2] = F[i2] - (dV_dpsi1 * dpsi1_dr2 + dV_dpsi2 * dpsi2_dr2 + dV_dr * r12_u);
  F[i2n] = F[i2n] - dV_dpsi2 * dpsi2_dr2n;
}

void eval_forces::fast_version_t::reset() {
  slices.clear();
  vect_slices.clear();
  fast_vect_slices.clear();
}

void eval_forces::fast_version_t::divide(
    std::vector<const simul::task *> &subtasks, int size_hint) {
  static constexpr auto V = eval_forces::elems_per_vect;
  auto fast_n = *super->fast_iter_end, fast_vect_n = fast_n / V;
  auto n = super->size(), vect_n = n / V;
  int vect_size_hint = size_hint / V;

  for (int fast_from = 0; fast_from < fast_vect_n;
       fast_from += vect_size_hint) {
    int fast_to = std::min(fast_from + vect_size_hint, fast_vect_n);
    fast_vect_slices.emplace_back(super, fast_from, fast_to);
  }

  for (int vect_from = fast_vect_n; vect_from < vect_n;
       vect_from += vect_size_hint) {
    int vect_to = std::min(vect_from + vect_size_hint, vect_n);
    vect_slices.emplace_back(super, vect_from, vect_to);
  }

  for (int from = vect_n * V; from < n; from += size_hint) {
    int to = std::min(from + size_hint, n);
    slices.emplace_back(super, from, to);
  }

  for (auto &slice_ : fast_vect_slices)
    subtasks.push_back(&slice_);

  for (auto &slice_ : vect_slices)
    subtasks.push_back(&slice_);

  for (auto &slice_ : slices)
    subtasks.push_back(&slice_);
}

eval_forces::fast_version_t::fast_version_t(const cg::pid::eval_forces *super)
    : super{super} {}

eval_forces::fast_version_t::slice::slice(const cg::pid::eval_forces *eval,
                                          int from, int to)
    : eval{eval}, from{from}, to{to} {}

void eval_forces::fast_version_t::slice::run() const {
  for (int idx = from; idx < to; ++idx)
    eval->iter(idx);
}

eval_forces::fast_version_t::fast_vect_slice::fast_vect_slice(
    const cg::pid::eval_forces *eval, int from, int to)
    : eval{eval}, from{from}, to{to} {}

void eval_forces::fast_version_t::fast_vect_slice::run() const {
  for (int idx = from; idx < to; ++idx)
    eval->fast_vect_iter(idx);
}

eval_forces::fast_version_t::vect_slice::vect_slice(
    const cg::pid::eval_forces *eval, int from, int to)
    : eval{eval}, from{from}, to{to} {}

void eval_forces::fast_version_t::vect_slice::run() const {
  for (int idx = from; idx < to; ++idx)
    eval->vect_iter(idx);
}

eval_forces::fast_version_t eval_forces::fast_version() const {
  return fast_version_t(this);
}

} // namespace cg::pid