#pragma once
#include "../nat_dih.h"

namespace cg::cnd {
class eval_forces {
public:
  real CDA, CDB;
  nitro::vector<vec3r> const *r;
  nitro::vector<vec3r> *F;
  nitro::vector<nat_dih> const *dihedrals;
  real *V;

public:
  template <typename E> void iter(nat_dih_expr<E> const &nat_dih) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::cnd
