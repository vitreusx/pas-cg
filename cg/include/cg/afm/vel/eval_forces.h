#pragma once
#include "tip.h"
#include <cg/base_forces/harmonic.h>

namespace cg::afm::vel {
class eval_forces {
public:
  harmonic afm_force;

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  real *t, *V;
  vect::view<tip> afm_tips;

public:
  template <typename E> void iter(tip_expr<E> &tip) const;
  void operator()() const;
  void omp_async() const;
  real compute_force(vel::tip const &tip) const;
};
} // namespace cg::afm::vel