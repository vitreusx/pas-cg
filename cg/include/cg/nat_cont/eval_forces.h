#pragma once
#include "nat_cont.h"
#include <cg/base_forces/disulfide.h>
#include <cg/types/box.h>

namespace cg::nat_cont {
class eval_forces {
public:
  real depth;
  std::optional<disulfide_force> disulfide;

public:
  nitro::const_view<vec3r> r;
  nitro::view<vec3r> F;
  box<real> const *simul_box;
  nitro::vector<nat_cont> const *contacts;
  real *V;

public:
  template <typename E> void iter(nat_cont_expr<E> const &nat_cont) const;
  void operator()() const;
  void omp_async() const;
};
} // namespace cg::nat_cont