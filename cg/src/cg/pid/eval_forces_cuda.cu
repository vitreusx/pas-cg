#include <cg/pid/eval_forces_cuda.h>

namespace cg::pid {
static __global__ void kernel(eval_forces_cuda eval) {
  auto idx = blockIdx.x * blockDim.x + blockIdx.x;
  if (idx >= eval.bundles.size())
    return;

  auto bundle = eval.bundles[idx];

  auto i1 = bundle.i1(), i2 = bundle.i2();
  if ((aa_code)eval.atype[i1] == aa_code::PRO ||
      (aa_code)eval.atype[i2] == aa_code::PRO)
    return;

  vec3r r1 = eval.r[i1], r2 = eval.r[i2];
  auto r12 = eval.simul_box.wrap<vec3r>(r1, r2);
  auto r12_rn = norm_inv(r12);
  if ((real)1.0 > r12_rn * eval.cutoff)
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
  auto n = eval.r.size();

  lambda const *bb_lam_1, *bb_lam_2;
  sink_lj const *bb_lj;

  {
    auto i1_ = i1, i2_ = i1n, i3_ = i1p, i4_ = i2;
    //    auto r1_ = r[i1_], r2_ = r[i2_], r3_ = r[i3_], r4_ = r[i4_];
    //    auto rij = r1_ - r2_, rkj = r3_ - r2_, rkl = r3_ - r4_;
    auto rij = eval.r[i1_] - (i2_ < n ? (vec3r)eval.r[i2_] : vec3r());
    auto rkj = i3_ >= 0 ? eval.r[i3_] - (i2_ < n ? (vec3r)eval.r[i2_] : vec3r())
                        : vec3r();
    auto rkl = (i3_ >= 0 ? (vec3r)eval.r[i3_] : vec3r()) - eval.r[i4_];

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
      bb_lam_1 = &eval.bb_plus_lam;
      bb_lj = &eval.bb_plus_lj;
    } else {
      bb_lam_1 = &eval.bb_minus_lam;
      bb_lj = &eval.bb_minus_lj;
    }
    if (!bb_lam_1->supp(psi1) && (i2 - i1 != 3 || m != 1)) {
      bb_lam_1 = &eval.bb_minus_lam;
      bb_lj = &eval.bb_minus_lj;
    }
  }

  {
    auto i1_ = i2, i2_ = i2n, i3_ = i2p, i4_ = i1;
    //    auto r1_ = r[i1_], r2_ = r[i2_], r3_ = r[i3_], r4_ = r[i4_];
    //    auto rij = r1_ - r2_, rkj = r3_ - r2_, rkl = r3_ - r4_;
    auto rij = eval.r[i1_] - (i2_ < n ? (vec3r)eval.r[i2_] : vec3r());
    auto rkj = i3_ >= 0 ? eval.r[i3_] - (i2_ < n ? (vec3r)eval.r[i2_] : vec3r())
                        : vec3r();
    auto rkl = (i3_ >= 0 ? (vec3r)eval.r[i3_] : vec3r()) - eval.r[i4_];

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
      bb_lam_2 = &eval.bb_plus_lam;
      bb_lj = &eval.bb_plus_lj;
    } else {
      bb_lam_2 = &eval.bb_minus_lam;
      bb_lj = &eval.bb_minus_lj;
    }
    if (!bb_lam_2->supp(psi2) && (i2 - i1 != 3 || m != 1)) {
      bb_lam_2 = &eval.bb_minus_lam;
      bb_lj = &eval.bb_minus_lj;
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

  sink_lj ss_sink_lj = eval.ss_ljs[type];
  if (eval.ss_lam.supp(psi1) && eval.ss_lam.supp(psi2)) {
    auto [lam1, deriv1] = eval.ss_lam(psi1);
    auto [lam2, deriv2] = eval.ss_lam(psi2);

    if (lam1 * lam2 > (real)5e-5 && ss_sink_lj.depth() > (real)5e-5) {
      auto [lj_V, lj_dV_dr] = ss_sink_lj(r12_n, r12_rn);

      V_ += lam1 * lam2 * lj_V;
      dV_dpsi1 += deriv1 * lam2 * lj_V;
      dV_dpsi2 += deriv2 * lam1 * lj_V;
      dV_dr += lam1 * lam2 * lj_dV_dr;
    }
  }

  atomicAdd(eval.V, V_);

  auto r12_u = r12 * r12_rn;
  if (i1p >= 0)
    eval.F[i1p].cudaAtomicSub(dV_dpsi1 * dpsi1_dr1p);
  eval.F[i1].cudaAtomicSub(dV_dpsi1 * dpsi1_dr1 + dV_dpsi2 * dpsi2_dr1 -
                           dV_dr * r12_u);
  if (i1n < n)
    eval.F[i1n].cudaAtomicSub(dV_dpsi1 * dpsi1_dr1n);
  if (i2p >= 0)
    eval.F[i2p].cudaAtomicSub(dV_dpsi2 * dpsi2_dr2p);
  eval.F[i2].cudaAtomicSub(dV_dpsi1 * dpsi1_dr2 + dV_dpsi2 * dpsi2_dr2 +
                           dV_dr * r12_u);
  if (i2n < n)
    eval.F[i2n].cudaAtomicSub(dV_dpsi2 * dpsi2_dr2n);
}

void eval_forces_cuda::operator()(int block_size, cudaStream_t stream) const {
  auto n = bundles.size();
  dim3 block(block_size), grid(n / block_size + 1);
  kernel<<<grid, block, 0, stream>>>(*this);
}
} // namespace cg::pid