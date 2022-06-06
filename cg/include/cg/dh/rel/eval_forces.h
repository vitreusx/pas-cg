#pragma once
#include "../pair.h"
#include <cg/types/box.h>

namespace cg::rel_dh {
class eval_forces {
public:
  real screen_dist_inv;
  real V_factor;
  void set_V_factor(real factor);

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  box<real> const *simul_box;
  vect::vector<dh::pair> const *es_pairs;
  real *V;

public:
  template <typename E> void iter(dh::pair_expr<E> const &pair) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::rel_dh