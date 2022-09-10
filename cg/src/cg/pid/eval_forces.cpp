#include <cg/pid/eval_forces.h>
#include <iostream>
#include <vcl/vectormath_trig.h>

namespace cg::pid {

// static auto acos_apx1(Vec4d x) {
//   static constexpr auto a = -0.939115566365855, b = 0.9217841528914573,
//                         c = -1.2845906244690837, d = 0.295624144969963174;
//   auto x2 = x * x;
//   return M_PI_2 + x * (a + b * x2) / (1.0 + x2 * (c + d * x2));
// }

void eval_forces::operator()() const {
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

  (void)0;
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

  auto r12 = norm(simul_box->wrap<vec3r>(r1, r2));
  auto r_min = ss_ljs[bundle.type()].r_high();
  auto sigma = r_min * C216_INV;
  return r12 <= active_thr * sigma;
}

void eval_forces::for_slice(int from, int to) const {
#if USE_VECTORIZED_IMPLS
#ifdef __AVX512F__
  static constexpr std::size_t W = 512;
#else
  static constexpr std::size_t W = 256;
#endif
  static constexpr std::size_t N = W / (8 * sizeof(real));

  int idx = from;
  for (; idx % N != 0 && idx < to; ++idx)
    iter(bundles->at(idx));
  for (; idx < to - (int)N + 1; idx += N)
    vect_iter<N, W>(idx / N);
  for (; idx < to; ++idx)
    iter(bundles->at(idx));
#else
  for (int idx = from; idx < to; ++idx)
    iter(bundles->at(idx));
#endif
}

int eval_forces::total_size() const {
  return bundles->size();
}

template <typename Mask, typename E, typename = void>
struct select_impl {
  static auto impl(Mask const &mask, E const &if_true, E const &if_false) {
    return ::select(mask, if_true, if_false);
  }
};

template <typename Mask, typename E>
auto select_(Mask const &mask, E const &if_true, E const &if_false) {
  return select_impl<Mask, E>::impl(mask, if_true, if_false);
}

template <typename Mask, typename E>
struct select_impl<Mask, E, std::enable_if_t<vect::is_indexed_v<E>>> {
  template <std::size_t... Idxes>
  static auto aux(Mask const &mask, E const &if_true, E const &if_false,
                  vect::ind_seq<Idxes...>) {
    return E(select_(mask, if_true.template get<Idxes>(),
                     if_false.template get<Idxes>())...);
  }

  static auto impl(Mask const &mask, E const &if_true, E const &if_false) {
    return aux(mask, if_true, if_false, vect::idxes_t<E>{});
  }
};

template <typename E>
auto norm_(vec3_expr<E> const &v) {
  return ::sqrt(norm_squared(v));
}

// template <typename E>
// std::ostream &operator<<(std::ostream &os, vec3_expr<E> const &v) {
//   os << "(" << v.x() << " " << v.y() << " " << v.z() << ")";
//   return os;
// }
//
// std::ostream &operator<<(std::ostream &os, vect::lane<vec3r, 4, 256> const
// &v) {
//   vect::vector<vec3r> V(4);
//   V.template at_lane<4, 256>(0) = v;
//   os << "[";
//   for (std::size_t idx = 0; idx < 4; ++idx) {
//     if (idx > 0)
//       os << " ";
//     os << V[idx];
//   }
//   os << "]";
//   return os;
// }
//
// template <typename T, std::size_t N>
// std::ostream &operator<<(std::ostream &os, std::array<T, N> const &A) {
//   os << "[";
//   for (std::size_t idx = 0; idx < N; ++idx) {
//     if (idx > 0)
//       os << " ";
//     os << A[idx];
//   }
//   os << "]";
//   return os;
// }

// std::ostream &operator<<(std::ostream &os, Vec4q v) {
//   std::array<int64_t, 4> A;
//   nitro::def::store(v, A.data());
//   return (os << A);
// }
//
// std::ostream &operator<<(std::ostream &os, Vec4uq v) {
//   std::array<uint64_t, 4> A;
//   nitro::def::store(v, A.data());
//   return (os << A);
// }
//
// std::ostream &operator<<(std::ostream &os, Vec4d v) {
//   std::array<double, 4> A;
//   nitro::def::store(v, A.data());
//   return (os << A);
// }
//
// std::ostream &operator<<(std::ostream &os, Vec4qb v) {
//   std::array<bool, 4> A;
//   nitro::def::store(v, A.data());
//   return (os << A);
// }
//
template <typename T, std::size_t N, std::size_t W>
auto to_array(vect::lane<T, N, W> const &lane) {
  vect::array<T, N> elems;
  elems.template at_lane<N, W>(0) = lane;
  return elems;
}

// #define PPRINT(var) std::cout << #var " = " << (var) << '\n';
// #define PPRINT(var) (void)0

template <std::size_t N, std::size_t W>
void eval_forces::vect_iter(int lane_idx) const {
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

} // namespace cg::pid