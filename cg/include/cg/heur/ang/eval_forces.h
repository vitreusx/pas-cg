#pragma once
#include "heur_angle.h"

namespace cg::heur_ang {
class eval_forces {
public:
  static constexpr int POLY_DEG = 6, NUM_TYPES = 9;
  real poly_coeffs[POLY_DEG + 1][NUM_TYPES];

public:
  nitro::const_view<vec3r> r;
  nitro::view<vec3r> F;
  nitro::const_view<heur_ang> angles;
  real *V;

public:
  template <typename E> void iter(heur_ang_expr<E> const &heur_ang) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::heur_ang