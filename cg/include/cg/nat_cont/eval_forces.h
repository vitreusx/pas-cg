#pragma once
#include "nat_cont.h"
#include <cg/base_forces/disulfide.h>
#include <cg/types/box.h>

namespace cg::nat_cont {
class eval_forces {
public:
  real depth, active_thr;
  std::optional<disulfide_force> disulfide;

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  box<real> const *simul_box;
  vect::vector<nat_cont> *contacts;
  vect::view<nat_cont> all_contacts;
  real *V, *t;

public:
  template <typename E> void iter(nat_cont_expr<E> &nat_cont) const;
  void operator()() const;
  void omp_async() const;

  bool is_active(nat_cont const &nc) const;
};
} // namespace cg::nat_cont