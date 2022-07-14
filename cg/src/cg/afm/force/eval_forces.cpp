#include <cg/afm/force/const.h>
namespace cg::afm::force {

void eval_forces::operator()() const {
  for (int idx = 0; idx < afm_tips.size(); ++idx) {
    iter(afm_tips[idx]);
  }
}

template <typename E> void eval_forces::iter(tip_expr<E> const &tip) const {
  F[tip.res_idx()] += tip.pull_force();
}

void eval_forces::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < afm_tips.size(); ++idx) {
    iter(afm_tips[idx]);
  }
}
} // namespace cg::afm::force