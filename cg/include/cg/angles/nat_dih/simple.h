#pragma once
#include "nat_dih.h"
#include <cg/simul/runtime.h>

namespace cg::snd {
class eval_forces : public simul::sliceable_task {
public:
  real CDH;
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  vect::const_view<nat_dih> dihedrals;
  real *V;

public:
  template <typename E> void iter(nat_dih_expr<E> const &nat_dih) const;
  void operator()() const;
  void omp_async() const;

  void for_slice(int from, int to) const override;
  int total_size() const override;

};
} // namespace cg::snd
