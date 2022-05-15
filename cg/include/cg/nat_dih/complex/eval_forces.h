#pragma once
#include "../nat_dih.h"

namespace cg::cnd {
class eval_forces {
public:
  real CDA, CDB;
  nitro::const_view<vec3r> r;
  nitro::view<vec3r> F;
  nitro::const_view<nat_dih> dihedrals;
  real *V;

public:
  template <typename E> void iter(nat_dih_expr<E> const &nat_dih) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::cnd