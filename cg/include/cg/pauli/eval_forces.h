#pragma once
#include "pair.h"
#include <cg/nl/data.h>
#include <cg/simul/runtime.h>

namespace cg::pauli {

class eval_forces : public simul::sliceable_task {
public:
  real depth, r_excl;

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  sbox::pbc<real> const *simul_box;
  vect::vector<pair> const *pairs;
  real *V;

public:
  template <typename E> void iter(pair_expr<E> const &pair) const;
  void operator()() const;
  void omp_async() const;

  void for_slice(int from, int to) const override;
  int total_size() const override;
  int slice_size() const override;
};
} // namespace cg::pauli