#pragma once
#include "bundle.h"
#include "lambda.h"
#include <cg/base_forces/lj.h>
#include <cg/base_forces/sink_lj.h>
#include <cg/sbox/pbc.h>

namespace cg::pid {
class eval_forces_cuda {
public:
  lambda bb_plus_lam, bb_minus_lam, ss_lam;
  sink_lj bb_plus_lj, bb_minus_lj;
  vect::const_view<sink_lj> ss_ljs;
  real cutoff;

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  sbox::pbc<real> simul_box;
  vect::const_view<bundle> bundles;
  real *V;
  vect::const_view<int> prev, next;
  vect::const_view<amino_acid> atype;

public:
  void operator()(int block_size, cudaStream_t stream = 0) const;
};
} // namespace cg::pid