#pragma once
#include "pair.h"
#include <cg/sbox/pbc.h>

namespace cg::const_dh {
class eval_forces_cuda {
public:
  real screen_dist_inv, V_factor, cutoff;
  void set_V_factor(real permittivity);

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  sbox::pbc<real> simul_box;
  vect::const_view<dh::pair> es_pairs;
  real *V;

public:
  void operator()(int block_size, cudaStream_t stream = 0) const;
};
} // namespace cg::const_dh