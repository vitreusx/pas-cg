#include "features/tether/eval_forces.h"
#include "base_forces/harmonic.h"
using namespace cg::tether;

void eval_forces::operator()() const {
  for (int idx = 0; idx < tethers->size(); ++idx) {
    iter(tethers->at(idx));
  }
}

void eval_forces::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < tethers->size(); ++idx) {
    iter(tethers->at(idx));
  }
}

template <typename E> void eval_forces::iter(pair_expr<E> const &tether) const {
  auto i1 = tether.i1(), i2 = tether.i2();
  auto nat_dist = tether.nat_dist();

  auto r1 = r->at(i1), r2 = r->at(i2);
  auto r12 = r2 - r1;

  auto r12_n = norm(r12);
  auto [V_, dV_dr] = harmonic(H1, H2, nat_dist)(r12_n);

  auto r12_u = r12 / r12_n;
  *V += V_;
  F->at(i1) += dV_dr * r12_u;
  F->at(i2) -= dV_dr * r12_u;
}