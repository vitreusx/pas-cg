#pragma once
#include "tip.h"

namespace cg::fafm {
class eval_forces {
public:
  nitro::vector<vec3r> *F;
  nitro::vector<tip> const *afm_tips;

public:
  template <typename E> void iter(tip_expr<E> const &tip) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::fafm