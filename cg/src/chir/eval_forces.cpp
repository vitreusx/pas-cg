#include "chir/eval_forces.h"
using namespace cg::chir;

void eval_forces::operator()() const {
  for (int idx = 0; idx < quads.size(); ++idx) {
    iter(quads[idx]);
  }
}

template <typename E>
void eval_forces::iter(chiral_quad_expr<E> const &quad) const {
  auto i1 = quad.i1(), i2 = quad.i2(), i3 = quad.i3(), i4 = quad.i4();
  auto nat_chir = quad.nat_chir();
  auto nat_factor = quad.nat_factor();

  auto r1 = r[i1], r2 = r[i2], r3 = r[i3], r4 = r[i4];
  auto r12 = r2 - r1, r23 = r3 - r2, r34 = r4 - r3;
  auto x12_23 = cross(r12, r23), x12_34 = cross(r12, r34),
       x23_34 = cross(r23, r34);

  auto chir = dot(r12, x23_34) * nat_factor;
  auto chir_diff = chir - nat_chir;

  *V += 0.5f * e_chi * chir_diff * chir_diff;

  auto f = e_chi * chir_diff * nat_factor;
  F[i1] += f * x12_23;
  F[i2] -= f * (x12_34 + x23_34);
  F[i3] += f * (x12_23 + x12_34);
  F[i4] -= f * x12_23;
}

void eval_forces::omp_async() const {
#pragma omp for schedule(static) nowait
  for (int idx = 0; idx < quads.size(); ++idx) {
    iter(quads[idx]);
  }
}