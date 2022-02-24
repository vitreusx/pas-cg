#pragma once
#include "heur_dih.h"

namespace cg::heur_dih {
class eval_forces {
public:
  static constexpr int NUM_TERMS = 6, NUM_TYPES = 9;
  real coeffs[NUM_TERMS][NUM_TYPES];

public:
  nitro::const_view<vec3r> r;
  nitro::view<vec3r> F;
  nitro::const_view<heur_dih> dihedrals;
  real *V;

public:
  template <typename E> void iter(heur_dih_expr<E> const &heur_dih) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::heur_dih