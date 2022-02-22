#pragma once
#include "tip.h"
#include <cg/base_forces/harmonic.h>

namespace cg::vafm {
class eval_forces {
public:
  harmonic afm_force;

public:
  nitro::vector<vec3r> const *r;
  nitro::vector<vec3r> *F;
  real *t;
  nitro::vector<tip> const *afm_tips;

public:
  template <typename E> void iter(tip_expr<E> const &tip) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::vafm