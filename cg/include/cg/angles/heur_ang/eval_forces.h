#pragma once
#include "heur_angle.h"
#include <cg/simul/runtime.h>

namespace cg::heur_ang {
class eval_forces : public simul::sliceable_task {
public:
  static constexpr int POLY_DEG = 6, NUM_TYPES = 9;
  real poly_coeffs[POLY_DEG + 1][NUM_TYPES];

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  vect::const_view<heur_ang> angles;
  real *V;

public:
  template <typename E> void iter(heur_ang_expr<E> const &heur_ang) const;
  void operator()() const;
  void omp_async() const;

  void for_slice(int from, int to) const override;
  int total_size() const override;

};
} // namespace cg::heur_ang