#include <cg/afm/force/eval_forces.h>
namespace cg::afm::force {

void eval_forces::iter(int idx) const {
  auto tip = afm_tips[idx];
  F[tip.res_idx()] += tip.pull_force();
  tip.avg_vel().add(*t, v[tip.res_idx()]);
}

int eval_forces::size() const {
  return afm_tips.size();
}

} // namespace cg::afm::force