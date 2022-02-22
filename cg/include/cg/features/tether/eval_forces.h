#pragma once
#include "pair.h"

namespace cg::tether {
class eval_forces {
public:
  real H1, H2;
  real def_length;

public:
  nitro::vector<vec3r> const *r;
  nitro::vector<vec3r> *F;
  nitro::vector<pair> const *tethers;
  real *V;

public:
  template <typename E> void iter(pair_expr<E> const &tether) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::tether