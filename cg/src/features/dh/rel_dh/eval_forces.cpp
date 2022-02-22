#include "features/dh/rel_dh/eval_forces.h"
using namespace cg::rel_dh;

void eval_forces::set_V_factor(real permittivity) {
  V_factor = 1.0 / (4.0 * M_PI * permittivity);
}

void eval_forces::operator()() const {
  for (int idx = 0; idx < es_pairs->size(); ++idx) {
    iter(es_pairs->at(idx));
  }
}

template <typename E> void eval_forces::iter(dh::pair_expr<E> const &es) const {
  auto i1 = es.i1(), i2 = es.i2();
  auto q1_x_q2 = es.q1_x_q2();

  auto r1 = r->at(i1), r2 = r->at(i2);
  auto r12 = box->r_uv(r1, r2);
  auto r12_n = norm(r12), r12_rn = 1.0f / r12_n;
  auto r12_u = r12 * r12_rn;

  auto Vij =
      V_factor * q1_x_q2 * exp(-r12_n * screen_dist_inv) * r12_rn * r12_rn;
  auto dVij_dr = -Vij * (screen_dist_inv + r12_rn);

  *V += Vij;

  auto f = r12_u * dVij_dr;
  F->at(i1) += f;
  F->at(i2) -= f;
}

void eval_forces::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < es_pairs->size(); ++idx) {
    iter(es_pairs->at(idx));
  }
}