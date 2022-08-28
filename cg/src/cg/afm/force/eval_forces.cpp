#include <cg/afm/force/eval_forces.h>
namespace cg::afm::force {

void eval_forces::operator()() const {
  for (int idx = 0; idx < afm_tips.size(); ++idx) {
    iter(afm_tips[idx]);
  }
}

template <typename E> void eval_forces::iter(tip_expr<E> &tip) const {
  F[tip.res_idx()] += tip.pull_force();
  tip.avg_vel().add(*t, v[tip.res_idx()]);
}

void eval_forces::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < afm_tips.size(); ++idx) {
    iter(afm_tips[idx]);
  }
}

void eval_forces::for_slice(int from, int to) const {
  for (int idx = from; idx < to; ++idx)
    iter(afm_tips[idx]);
}

int eval_forces::total_size() const {
  return afm_tips.size();
}


} // namespace cg::afm::force