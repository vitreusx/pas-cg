#pragma once
#include "heur_dih.h"

namespace cg::heur_dih {
class eval_forces {
public:
  struct coeffs_t {
    real const_[aa_heur_pair::NUM_TYPES];
    real sin[aa_heur_pair::NUM_TYPES];
    real cos[aa_heur_pair::NUM_TYPES];
    real sin2[aa_heur_pair::NUM_TYPES];
    real cos2[aa_heur_pair::NUM_TYPES];
    real sin_cos[aa_heur_pair::NUM_TYPES];
  };
  coeffs_t coeffs;

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