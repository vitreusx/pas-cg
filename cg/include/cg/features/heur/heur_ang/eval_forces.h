#pragma once
#include "heur_angle.h"

namespace cg::heur_ang {
class eval_forces {
public:
  static constexpr int POLY_DEG = 6, NUM_TYPES = 9;
  real poly_coeffs[POLY_DEG + 1][NUM_TYPES];

public:
  nitro::vector<vec3r> const *r;
  nitro::vector<vec3r> *F;
  nitro::vector<heur_ang> const *angles;
  real *V;

public:
  template <typename E> void iter(heur_ang_expr<E> const &heur_ang) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::heur_ang