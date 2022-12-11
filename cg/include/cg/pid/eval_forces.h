#pragma once
#include "bundle.h"
#include "lambda.h"
#include <cg/base_forces/lj.h>
#include <cg/base_forces/sink_lj.h>
#include <cg/sbox/pbc.h>
#include <cg/simul/runtime.h>
#include <cg/simul/sched.h>

namespace cg::pid {
class eval_forces : public simul::vect_iter_divisible_mixin<eval_forces> {
public:
  lambda bb_plus_lam, bb_minus_lam, ss_lam;
  sink_lj bb_plus_lj, bb_minus_lj;
  vect::vector<sink_lj> ss_ljs;
  real cutoff, active_thr;

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  sbox::pbc<real> const *simul_box;
  vect::vector<bundle> const *bundles;
  int const *fast_iter_end;
  real *V, *total_disp;
  vect::const_view<int> prev, next;
  vect::const_view<amino_acid> atype;

public:
  bool is_active(bundle const &bundle) const;

public:
  void iter(int idx) const;
  void vect_iter(int vect_idx) const;
  int size() const;

  static inline constexpr int elems_per_vect = vect::VECT_BYTES / sizeof(real);
};
} // namespace cg::pid