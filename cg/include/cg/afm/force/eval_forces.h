#pragma once
#include "tip.h"
#include <cg/simul/runtime.h>

namespace cg::afm::force {
class eval_forces : public simul::sliceable_task {
public:
  vect::view<vec3r> F;
  vect::view<tip> afm_tips;
  vect::const_view<vec3r> v;
  real const *t;

public:
  template <typename E> void iter(tip_expr<E> &tip) const;
  void operator()() const;
  void omp_async() const;

  void for_slice(int from, int to) const override;
  int total_size() const override;
  int slice_size() const override;
};
} // namespace cg::afm::force