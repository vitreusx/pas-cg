#include <cg/dh/rel.h>
namespace cg::rel_dh {

void eval_forces::set_V_factor(real permittivity) {
  V_factor = 1.0 / (4.0 * M_PI * permittivity);
}

void eval_forces::operator()() const {
  if (!es_pairs)
    return;
  
  for (int idx = 0; idx < es_pairs->size(); ++idx) {
    iter(es_pairs->at(idx));
  }
}

template <typename E>
void eval_forces::iter(dh::pair_expr<E> const &es) const {
  auto i1 = es.i1(), i2 = es.i2();
  auto q1_x_q2 = es.q1_x_q2();

  auto r1 = r[i1], r2 = r[i2];
  auto r12 = simul_box->wrap(r1, r2);

  auto r12_n = norm(r12);
  if (r12_n > cutoff)
    return;
  auto r12_rn = 1.0f / r12_n;
  auto r12_u = r12 * r12_rn;

  auto Vij =
      V_factor * q1_x_q2 * exp(-r12_n * screen_dist_inv) * r12_rn * r12_rn;
  auto dVij_dr = -Vij * (screen_dist_inv + (real)2.0 * r12_rn);

  *V += Vij;

  auto f = r12_u * dVij_dr;
  F[i1] += f;
  F[i2] -= f;
}

void eval_forces::omp_async() const {
  if (!es_pairs)
    return;

#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < es_pairs->size(); ++idx) {
    iter(es_pairs->at(idx));
  }
}

void eval_forces::for_slice(int from, int to) const {
  for (int idx = from; idx < to; ++idx)
    iter(es_pairs->at(idx));
}

int eval_forces::total_size() const {
  return es_pairs ? es_pairs->size() : 0;
}

} // namespace cg::rel_dh