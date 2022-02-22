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
  nitro::vector<sink_lj> const *ss_ljs;

public:
  nitro::vector<vec3r> const *r;
  nitro::vector<vec3r> *F;
  box<real> const *box;
  nitro::vector<bundle> const *bundles;
  real *V;
  nitro::vector<int> const *prev, *next;

public:
  template <typename E> void iter(bundle_expr<E> const &bundle) const;
  void operator()();
  void omp_async() const;
};
} // namespace cg::pid