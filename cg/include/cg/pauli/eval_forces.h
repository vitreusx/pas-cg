#pragma once
#include "pair.h"
#include <cg/nl/data.h>
namespace cg::pauli {

class eval_forces {
public:
  real depth, r_excl;

public:
  nitro::const_view<vec3r> r;
  nitro::view<vec3r> F;
  box<real> const *box;
  nitro::vector<pair> const *pairs;
  real *V;

public:
  template <typename E> void iter(pair_expr<E> const &pair) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::pauli