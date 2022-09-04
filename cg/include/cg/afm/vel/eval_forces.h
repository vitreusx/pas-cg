#pragma once
#include "tip.h"
#include <cg/base_forces/harmonic.h>
#include <cg/simul/runtime.h>

namespace cg::afm::vel {
class eval_forces : public simul::sliceable_task {
public:
  harmonic afm_force;

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  real *t, *V;
  vect::view<tip> afm_tips;

public:
  template <typename E> void iter(tip_expr<E> const &tip) const;
  void operator()() const;
  void omp_async() const;
  real compute_force(vel::tip const &tip) const;

  void for_slice(int from, int to) const override;
  int total_size() const override;
};
} // namespace cg::afm::vel