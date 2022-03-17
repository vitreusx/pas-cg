#include "qa/process_contacts.h"
namespace cg::qa {

real process_contacts::saturation_value(const contact &cont) const {
  auto saturation = min(*t - cont.ref_time(), cycle_time) * cycle_time_inv;
  if (cont.status() == BREAKING)
    saturation = 1.0f - saturation;
  return saturation;
}

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
  auto status = contact.status();

  auto r1 = r[i1], r2 = r[i2];
  auto r12 = simul_box->r_uv(r1, r2);
  auto r12_rn = norm_inv(r12);
  auto r12_u = r12 * r12_rn;

  auto saturation = saturation_value(contact);

  static contact_type ss_type =
      contact_type::SIDE_SIDE(aa_code::CYS, aa_code::CYS);
  if (!disulfide.has_value() || (int16_t)type != (int16_t)ss_type) {
    auto lj_force = ljs[type];
    auto r12_n = r12_rn * norm_squared(r12);
    auto [Vij, dVij_dr] = lj_force(r12_n, r12_rn);
    *V += saturation * Vij;
    auto f = saturation * dVij_dr * r12_u;
    F[i1] += f;
    F[i2] -= f;

    if (status == FORMING_OR_FORMED && saturation == 1.0f) {
      if (factor * lj_force.r_max() * r12_rn < 1.0f) {
        contact.status() = BREAKING;
        contact.ref_time() = *t;
      }
    } else if (status == BREAKING && saturation == 0.0f) {
#pragma omp critical
      removed->push_back(idx);
    }
  } else {
    auto [Vij, dVij_dr] = disulfide.value()(r12_rn);
    *V += saturation * Vij;
    auto f = saturation * dVij_dr * r12_u;
    F[i1] += f;
    F[i2] -= f;

    if (status == FORMING_OR_FORMED && saturation == 1.0f) {
      if (disulfide_special_criteria) {
        real r12_n = r12_rn * norm_squared(r12);
        if (abs(r12_n - ss_def_dist) > ss_dist_max_div &&
            neigh[i1] + neigh[i2] < max_neigh_count) {
          contact.status() = BREAKING;
          contact.ref_time() = *t;
        }
      } else {
        if (factor * ljs[ss_type].r_max() * r12_rn < 1.0f) {
          contact.status() = BREAKING;
          contact.ref_time() = *t;
        }
      }
    } else if (status == BREAKING && saturation == 0.0f) {
#pragma omp critical
      removed->push_back(idx);
    }
  }
}

void process_contacts::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < contacts->size(); ++idx) {
    iter(idx);
  }
}

void process_contacts::set_factor(double breaking_factor) {
  factor = breaking_factor * pow(2.0, -1.0 / 6.0);
}
} // namespace cg::qa