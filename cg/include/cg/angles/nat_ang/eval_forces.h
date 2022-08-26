#pragma once
#include "nat_ang.h"
#include <cg/simul/runtime.h>

namespace cg::nat_ang {
class eval_forces : public simul::sliceable_task {
public:
  real CBA;

public:
  real *V;
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  vect::const_view<nat_ang> angles;

public:
  template <typename E> void iter(nat_ang_expr<E> const &angle) const;
  void operator()() const;
  void omp_async() const;

  void for_slice(int from, int to) const override;
  int total_size() const override;
  int slice_size() const override;
};
} // namespace cg::nat_ang