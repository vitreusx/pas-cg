#include "nat_cont/eval_forces.h"
#include "base_forces/lj.h"
using namespace cg::nat_cont;

void eval_forces::operator()() const {
  for (int idx = 0; idx < contacts->size(); ++idx) {
    iter(contacts->at(idx));
  }
}

template <typename E>
void eval_forces::iter(nat_cont_expr<E> const &nat_cont) const {
  auto i1 = nat_cont.i1(), i2 = nat_cont.i2();
  auto nat_dist = nat_cont.nat_dist();

  auto r1 = r[i1], r2 = r[i2];
  auto r12 = box->r_uv(r1, r2);
  auto r12_rn = norm_inv(r12);

  auto r12_u = r12 * r12_rn;
  auto [V_, dV_dr] = lj(depth, nat_dist)(r12_rn);

  *V += V_;
  F[i1] += r12_u * dV_dr;
  F[i2] -= r12_u * dV_dr;
}

void eval_forces::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < contacts->size(); ++idx) {
    iter(contacts->at(idx));
  }
}