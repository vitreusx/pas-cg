#pragma once
#include "tip.h"

namespace cg::afm::force {
class eval_forces {
public:
  vect::view<vec3r> F;
  vect::view<tip> afm_tips;
  vect::const_view<vec3r> v;
  real const *t;

public:
  template <typename E> void iter(tip_expr<E> &tip) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::afm::force