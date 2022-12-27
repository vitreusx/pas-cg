#include <cg/dh/const.h>

namespace cg::const_dh {

void eval_forces::set_V_factor(real factor) {
  V_factor = 1.0 / (4.0 * M_PI * factor);
}

void eval_forces::iter(int idx) const {
  auto es = (*es_pairs)[idx];
  auto i1 = es.i1(), i2 = es.i2();
  auto q1_x_q2 = es.q1_x_q2();

  auto r1 = r[i1], r2 = r[i2];
  auto r12 = simul_box->wrap<vec3r>(r1, r2);

  auto r12_n = norm(r12);
  if (r12_n > cutoff)
    return;
  auto r12_rn = 1.0f / r12_n;
  auto r12_u = r12 * r12_rn;

  auto Vij = V_factor * q1_x_q2 * exp(-r12_n * screen_dist_inv) * r12_rn;
  auto dVij_dr = -Vij * (screen_dist_inv + r12_rn);

  *V += Vij;

  auto f = r12_u * dVij_dr;
  F[i1] += f;
  F[i2] -= f;
}

template <typename T, std::size_t N, std::size_t W>
static auto to_array(vect::lane<T, N, W> const &lane) {
  vect::array<T, N> elems;
  elems.template at_lane<N, W>(0) = lane;
  return elems;
}

void eval_forces::vect_iter(int idx) const {
  static constexpr auto N = elems_per_vect, W = vect::VECT_BITS;

  auto es = es_pairs->at_lane<N, W>(idx);
  auto i1 = es.i1(), i2 = es.i2();
  auto q1_x_q2 = es.q1_x_q2();

  auto r1 = r[i1], r2 = r[i2];
  auto r12 = simul_box->wrap<vect::lane<vec3r, N, W>>(r1, r2);

  auto r12_n = norm(r12);
  auto mask = r12_n <= cutoff;
  auto zero_v = vect::lane<real, N, W>((real)0.0);
  q1_x_q2 = select(mask, q1_x_q2, zero_v);

  auto r12_rn = (real)1.0 / r12_n;
  auto r12_u = r12 * r12_rn;

  auto Vij = V_factor * q1_x_q2 * exp(-r12_n * screen_dist_inv) * r12_rn;
  auto dVij_dr = -Vij * (screen_dist_inv + r12_rn);

  *V += horizontal_add(Vij);

  auto f = r12_u * dVij_dr;

  auto f_ = to_array<vec3r, N, W>(f);
  auto i1_ = to_array<int, N, W>(i1), i2_ = to_array<int, N, W>(i2);
  for (int x = 0; x < N; ++x) {
    auto f_x = f_[x];
    F[i1_[x]] += f_x;
    F[i2_[x]] -= f_x;
  }
}

int eval_forces::size() const {
  return es_pairs ? es_pairs->size() : 0;
}

} // namespace cg::const_dh