#include "features/force_afm/eval_forces.h"
using namespace cg::fafm;

void eval_forces::operator()() const {
  for (int idx = 0; idx < afm_tips->size(); ++idx) {
    iter(afm_tips->at(idx));
  }
}

template <typename E> void eval_forces::iter(tip_expr<E> const &tip) const {
  F->at(tip.res_idx()) += tip.pull_force();
}

void eval_forces::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < afm_tips->size(); ++idx) {
    iter(afm_tips->at(idx));
  }
}