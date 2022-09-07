#pragma once
#include "bundle.h"
#include "lambda.h"
#include <cg/base_forces/lj.h>
#include <cg/base_forces/sink_lj.h>
#include <cg/sbox/pbc.h>
#include <cg/simul/runtime.h>

namespace cg::pid {
class eval_forces : public simul::sliceable_task {
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
  real *V, *total_disp;
  vect::const_view<int> prev, next;
  vect::const_view<amino_acid> atype;

public:
  template <typename E>
  void iter(bundle_expr<E> const &bundle) const;
  void operator()() const;
  void omp_async() const;

  template<std::size_t N, std::size_t W>
  void vect_iter(int lane_idx) const;

  bool is_active(bundle const &bundle) const;

  void for_slice(int from, int to) const override;
  int total_size() const override;
};
} // namespace cg::pid