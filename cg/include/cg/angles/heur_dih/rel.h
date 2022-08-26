#pragma once
#include "heur_dih.h"
#include <cg/simul/runtime.h>

namespace cg::heur_dih {
class eval_forces : public simul::sliceable_task {
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
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  vect::const_view<heur_dih> dihedrals;
  real *V;

public:
  template <typename E> void iter(heur_dih_expr<E> const &heur_dih) const;
  void operator()() const;
  void omp_async() const;

  void for_slice(int from, int to) const override;
  int total_size() const override;
  int slice_size() const override;
};
} // namespace cg::heur_dih