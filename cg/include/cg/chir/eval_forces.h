#pragma once
#include "chiral_quad.h"
#include <cg/simul/runtime.h>
#include <cg/vect/vect.h>

namespace cg::chir {
class eval_forces : public simul::sliceable_task {
public:
  real e_chi;

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  vect::const_view<chiral_quad> quads;
  real *V;

public:
  template <typename E> void iter(chiral_quad_expr<E> const &quad) const;
  void operator()() const;
  void omp_async() const;

  void for_slice(int from, int to) const override;
  int total_size() const override;
  int slice_size() const override;
};
} // namespace cg::chir