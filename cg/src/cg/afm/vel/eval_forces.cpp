#include <cg/afm/vel/eval_forces.h>

namespace cg::afm::vel {
void eval_forces::iter(int idx) const {
  auto tip = afm_tips[idx];
  auto r_ = r[tip.res_idx()];
  auto cur_afm_pos = tip.afm_orig() + (*t - tip.ref_t()) * tip.afm_vel();

  auto r_afm = cur_afm_pos - r_;
  auto r_afm_n = norm(r_afm);
  if (r_afm_n > (real)1e-10) {
    auto r_afm_u = r_afm / r_afm_n;
    auto [V_, dV_dr] = afm_force(r_afm_n);
    *V += V_;
    auto force = dV_dr * r_afm_u;
    F[tip.res_idx()] += force;

    tip.avg_force().add(*t, -force);

    auto v = unit(tip.afm_vel());
    auto perp_force =
        force.x() * v.z() / (v.x() + v.z() * v.z() / v.x()) -
        force.z() / ((real)1.0 + (v.z() / v.x()) * (v.z() / v.x()));
    tip.avg_perp_force().add(*t, perp_force);
  }
}

int eval_forces::size() const {
  return afm_tips.size();
}

cg::real eval_forces::compute_force(const vel::tip &tip) const {
  auto r_ = r[tip.res_idx()];
  auto cur_afm_pos = tip.afm_orig() + (*t - tip.ref_t()) * tip.afm_vel();
  auto r_afm = cur_afm_pos - r_;
  auto r_afm_n = norm(r_afm);
  auto [_, dV_dr] = afm_force(r_afm_n);
  return dV_dr;
}
} // namespace cg::afm::vel