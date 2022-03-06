#include "afm/vel/eval_forces.h"
using namespace cg::afm::vel;

void eval_forces::operator()() const {
  for (int idx = 0; idx < afm_tips.size(); ++idx) {
    iter(afm_tips[idx]);
  }
}

template <typename E> void eval_forces::iter(tip_expr<E> const &tip) const {
  auto r_ = r[tip.res_idx()];
  auto cur_afm_pos = tip.afm_orig() + *t * tip.afm_vel();

  auto r_afm = cur_afm_pos - r_;
  auto r_afm_n = norm(r_afm), r_afm_rn = 1.0f / r_afm_n;
  auto r_afm_u = r_afm * r_afm_rn;

  auto [_, dV_dr] = afm_force(r_afm_n);
  F[tip.res_idx()] += dV_dr * r_afm_u;
}

void eval_forces::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < afm_tips.size(); ++idx) {
    iter(afm_tips[idx]);
  }
}