#pragma once
#include "pair.h"
#include <cg/simul/runtime.h>

namespace cg::tether {
class eval_forces : public simul::sliceable_task {
public:
  real H1, H2;
  real def_length;

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  vect::const_view<pair> tethers;
  real *V;

public:
  template <typename E> void iter(pair_expr<E> const &tether) const;
  void operator()() const;
  void omp_async() const;

  void for_slice(int from, int to) const override;
  int total_size() const override;

};
} // namespace cg::tether