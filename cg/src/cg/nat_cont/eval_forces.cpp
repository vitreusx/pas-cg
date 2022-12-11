#include <cg/base_forces/lj.h>
#include <cg/nat_cont/eval_forces.h>

namespace cg::nat_cont {
void eval_forces::iter(int idx) const {
  auto nat_cont = contacts->at(idx);
  auto i1 = nat_cont.i1(), i2 = nat_cont.i2();
  auto nat_dist = nat_cont.nat_dist();

  auto r1 = r[i1], r2 = r[i2];
  auto r12 = simul_box->wrap<vec3r>(r1, r2);
  auto r12_n = norm(r12);

  if (r12_n >= cutoff)
    return;

  auto r12_u = r12 / r12_n;

  real V_, dV_dr;
  if (disulfide.has_value() && nat_cont.is_ssbond()) {
    std::tie(V_, dV_dr) = disulfide.value()(r12_n);
  } else {
    auto r12_rn = (real)1.0 / r12_n;
    std::tie(V_, dV_dr) = lj(depth, nat_dist)(r12_rn);
  }
  dV_dr = clamp(dV_dr, (real)-1e3, (real)1e3);

  *V += V_;
  F[i1] += r12_u * dV_dr;
  F[i2] -= r12_u * dV_dr;

  bool active_now = r12_n <= nat_dist * C216_INV * breaking_threshold;
  if (active_now ^ nat_cont.active()) {
    auto ref_idx = nat_cont.all_cont_idx();
    decltype(auto) ref_contact = all_contacts[ref_idx];
    nat_cont.active() = ref_contact.active() = active_now;

    if (nat_cont.change_t() < 0) {
      nat_cont.change_t() = ref_contact.change_t() = *t;
#pragma omp critical
      ++*num_changed;
    }
  }
}

int eval_forces::size() const {
  return contacts->size();
}

} // namespace cg::nat_cont