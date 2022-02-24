#pragma once
#include "bundle.h"
#include "lambda.h"
#include <cg/base_forces/lj.h>
#include <cg/base_forces/sink_lj.h>
#include <cg/types/box.h>

namespace cg::pid {
class eval_forces {
public:
  lambda bb_plus_lam, bb_minus_lam, ss_lam;
  lj bb_plus_lj, bb_minus_lj;
  nitro::const_view<sink_lj> ss_ljs;

public:
  nitro::const_view<vec3r> r;
  nitro::view<vec3r> F;
  box<real> const *box;
  nitro::vector<bundle> const *bundles;
  real *V;
  nitro::const_view<int> prev, next;

public:
  template <typename E> void iter(bundle_expr<E> const &bundle) const;
  void operator()();
  void omp_async() const;
};
} // namespace cg::pid