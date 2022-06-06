#pragma once
#include "nat_ang.h"

namespace cg::nat_ang {
class eval_forces {
public:
  real k;

public:
  real *V;
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  vect::const_view<nat_ang> angles;

public:
  template <typename E> void iter(nat_ang_expr<E> const &angle) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::nat_ang