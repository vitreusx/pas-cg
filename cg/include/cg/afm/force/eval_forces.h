#pragma once
#include "tip.h"

namespace cg::afm::force {
class eval_forces {
public:
  nitro::view<vec3r> F;
  nitro::const_view<tip> afm_tips;

public:
  template <typename E> void iter(tip_expr<E> const &tip) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::afm::force