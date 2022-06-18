#pragma once
#include "pair.h"

namespace cg::local_rep {
class eval_forces {
public:
  real depth, cutoff;

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  vect::const_view<pair> pairs;
  real *V;

public:
  template <typename E> void iter(pair_expr<E> const &pair) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::local_rep