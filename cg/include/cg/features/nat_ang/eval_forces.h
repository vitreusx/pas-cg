#pragma once
#include "nat_ang.h"

namespace cg::nat_ang {
class eval_forces {
public:
  real k;

public:
  real *V;
  nitro::vector<vec3r> const *r;
  nitro::vector<vec3r> *F;
  nitro::vector<nat_ang> const *angles;

public:
  template <typename E> void iter(nat_ang_expr<E> const &angle) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::nat_ang