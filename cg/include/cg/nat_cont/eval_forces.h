#pragma once
#include "nat_cont.h"
#include <cg/base_forces/disulfide.h>
#include <cg/sbox/pbc.h>

namespace cg::nat_cont {
class eval_forces {
public:
  real depth, breaking_threshold, cutoff;
  std::optional<disulfide_force> disulfide;

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  sbox::pbc<real> const *simul_box;
  vect::vector<nat_cont> *contacts;
  vect::view<nat_cont> all_contacts;
  real *V, *t;

public:
  template <typename E> void iter(nat_cont_expr<E> &nat_cont) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::nat_cont