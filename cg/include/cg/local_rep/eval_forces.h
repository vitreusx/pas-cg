#pragma once
#include "pair.h"
#include <cg/simul/runtime.h>

namespace cg::local_rep {
class eval_forces : public simul::sliceable_task {
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

  void for_slice(int from, int to) const override;
  int total_size() const override;
  int slice_size() const override;
};
} // namespace cg::local_rep