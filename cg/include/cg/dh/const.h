#pragma once
#include "pair.h"
#include <cg/sbox/pbc.h>
#include <cg/simul/runtime.h>
#include <optional>

namespace cg::const_dh {
class eval_forces : public simul::sliceable_task {
public:
  real screen_dist_inv, V_factor, cutoff;
  void set_V_factor(real permittivity);

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  sbox::pbc<real> const *simul_box;
  vect::const_view<dh::pair> es_pairs;
  real *V;

public:
  template <typename E> void iter(dh::pair_expr<E> const &pair) const;
  void operator()() const;
  void omp_async() const;

  void for_slice(int from, int to) const override;
  int total_size() const override;

};
} // namespace cg::const_dh