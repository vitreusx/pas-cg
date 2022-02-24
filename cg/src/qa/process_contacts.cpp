#include "qa/process_contacts.h"
using namespace cg::qa;

void process_contacts::operator()() const {
  for (int idx = 0; idx < contacts->size(); ++idx) {
    iter(idx);
  }
}

void process_contacts::iter(int idx) const {
  if (contacts->is_vacant(idx))
    return;

  auto contact = contacts->at(idx);

  auto i1 = contact.i1(), i2 = contact.i2();
  auto type = contact.type();
  auto status = contact.status();
  auto ref_time = contact.ref_time();
  auto sync_diff1 = contact.sync_diff1();
  auto sync_diff2 = contact.sync_diff2();

  auto r1 = r[i1], r2 = r[i2];
  auto r12 = box->r_uv(r1, r2);
  auto r12_rn = norm_inv(r12);
  auto r12_u = r12 * r12_rn;

  auto saturation = min(*t - ref_time, cycle_time) * cycle_time_inv;
  if (status == BREAKING)
    saturation = 1.0f - saturation;

  auto lj_force = ljs[type];
  auto r12_n = r12_rn * norm_squared(r12);
  auto [Vij, dVij_dr] = lj_force(r12_n, r12_rn);
  *V += saturation * Vij;
  auto f = saturation * dVij_dr * r12_u;
  F[i1] += f;
  F[i2] -= f;

  if (status == FORMING_OR_FORMED && saturation == 1.0f) {
    if (factor * lj_force.r_min() * r12_rn < 1.0f) {
      contact.status() = BREAKING;
      contact.ref_time() = *t;
    }
  } else if (status == BREAKING && saturation == 0.0f) {
#pragma omp critical
    {
      contacts->remove(idx);
      sync[i1] += sync_diff1;
      sync[i2] += sync_diff2;
      free_pairs->emplace_back(i1, i2);
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