#include <cg/base_forces/lj.h>
#include <cg/nat_cont/eval_forces.h>
namespace cg::nat_cont {

void eval_forces::operator()() const {
  for (int idx = 0; idx < contacts->size(); ++idx) {
    iter(contacts->at(idx));
  }
}

template <typename E> void eval_forces::iter(nat_cont_expr<E> &nat_cont) const {
  auto i1 = nat_cont.i1(), i2 = nat_cont.i2();
  auto nat_dist = nat_cont.nat_dist();

  auto r1 = r[i1], r2 = r[i2];
  auto r12 = simul_box->wrap(r1, r2);
  auto r12_rn = norm_inv(r12);

  auto r12_u = r12 * r12_rn;

  real V_, dV_dr;
  if (disulfide.has_value() && nat_cont.is_ssbond()) {
    auto r12_n = (real)1 / r12_rn;
    std::tie(V_, dV_dr) = disulfide.value()(r12_n);
  } else {
    std::tie(V_, dV_dr) = lj(depth, nat_dist)(r12_rn);
  }

  *V += V_;
  F[i1] += r12_u * dV_dr;
  F[i2] -= r12_u * dV_dr;

  auto r12_n = (real)1.0 / r12_rn;
  if (r12_n <= nat_dist * active_thr) {
    if (!nat_cont.formed()) {
      auto ref_idx = nat_cont.all_cont_idx();
      auto ref_contact = all_contacts[ref_idx];
      nat_cont.formed() = ref_contact.formed() = true;
      nat_cont.formation_t() = ref_contact.formation_t() = *t;
    }
  }
}

void eval_forces::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < contacts->size(); ++idx) {
    iter(contacts->at(idx));
  }
}

bool eval_forces::is_active(const nat_cont &nc) const {
  auto i1 = nc.i1(), i2 = nc.i2();
  auto nat_dist = nc.nat_dist();
  auto r1 = r[i1], r2 = r[i2];
  auto r12 = simul_box->wrap(r1, r2);
  return norm(r12) <= nat_dist * active_thr;
}
} // namespace cg::nat_cont