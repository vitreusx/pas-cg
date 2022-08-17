#include <cg/qa/process_contacts.h>
#include <cg/utils/math.h>

namespace cg::qa {

void process_contacts::operator()() const {
  for (int idx = 0; idx < contacts->size(); ++idx) {
    iter(idx);
  }
}

void process_contacts::iter(int idx) const {
  auto const &node = contacts->at(idx);

  if (node.is_vacant())
    return;

  auto contact = node.item();

  auto i1 = contact.i1(), i2 = contact.i2();
  auto type = contact.type();
  auto saturation = contact.saturation();

  auto r1 = r[i1], r2 = r[i2];
  auto r12 = simul_box->wrap(r1, r2);
  auto r12_rn = norm_inv(r12);
  if (fixed_cutoff.has_value() && r12_rn * fixed_cutoff.value() < 1.0)
    return;

  auto r12_u = r12 * r12_rn;

  static contact_type ss_type =
      contact_type::SIDE_SIDE(aa_code::CYS, aa_code::CYS);

  bool breaking = false;
  if (!disulfide.has_value() || (int16_t)type != (int16_t)ss_type) {
    auto lj_force = ljs[(int16_t)type];
    auto r12_n = r12_rn * norm_squared(r12);
    auto [Vij, dVij_dr] = lj_force(r12_n, r12_rn);
    *V += saturation * Vij;
    auto f = saturation * dVij_dr * r12_u;
    F[i1] += f;
    F[i2] -= f;

    breaking = r12_n > factor * lj_force.r_high();
  } else {
    auto [Vij, dVij_dr] = disulfide.value()(r12_rn);
    *V += saturation * Vij;
    auto f = saturation * dVij_dr * r12_u;
    F[i1] += f;
    F[i2] -= f;

    if (disulfide_special_criteria) {
      real r12_n = r12_rn * norm_squared(r12);
      breaking = abs(r12_n - ss_def_dist) > ss_dist_max_div &&
                 neigh[i1] + neigh[i2] < max_neigh_count;
    } else {
      breaking = factor * ljs[(int16_t)ss_type].r_high() * r12_rn < (real)1.0;
    }
  }

  if (breaking) {
    saturation = min(saturation - saturation_speed * dt, (real)0.0);
    if (saturation == (real)0.0) {
#pragma omp critical
      removed->push_back(idx);
    }
  } else {
    saturation = max(saturation + saturation_speed * dt, (real)0.0);
  }
  contact.saturation() = saturation;
}

void process_contacts::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < contacts->size(); ++idx) {
    iter(idx);
  }
}

void process_contacts::set_factor(double breaking_factor) {
  factor = breaking_factor * C216_INV;
}
} // namespace cg::qa