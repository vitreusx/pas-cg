#pragma once
#include "chiral_quad.h"
#include "nitro/nitro.h"

namespace cg::chir {
class eval_forces {
public:
  real e_chi;

public:
  nitro::vector<vec3r> const *r;
  nitro::vector<vec3r> *F;
  nitro::vector<chiral_quad> const *quads;
  real *V;

public:
  template <typename E> void iter(chiral_quad_expr<E> const &quad) const;

  void operator()() const;
  void omp_async() const;
};
} // namespace cg::chir