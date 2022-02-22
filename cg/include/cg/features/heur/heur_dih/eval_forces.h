#pragma once
#include "heur_dih.h"

namespace cg::heur_dih {
class eval_forces {
public:
  static constexpr int NUM_TERMS = 6, NUM_TYPES = 9;
  real coeffs[NUM_TERMS][NUM_TYPES];

public:
  nitro::vector<vec3r> const *r;
  nitro::vector<vec3r> *F;
  nitro::vector<heur_dih> const *dihedrals;
  real *V;

public:
  template <typename E> void iter(heur_dih_expr<E> const &heur_dih) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::heur_dih