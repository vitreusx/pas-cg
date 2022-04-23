#pragma once
#include "chiral_quad.h"
#include <cg/vect/vect.h>

namespace cg::chir {
class eval_forces {
public:
  real e_chi;

public:
  nitro::const_view<vec3r> r;
  nitro::view<vec3r> F;
  nitro::const_view<chiral_quad> quads;
  real *V;

public:
  template <typename E> void iter(chiral_quad_expr<E> const &quad) const;

  void operator()() const;
  void omp_async() const;
};
} // namespace cg::chir