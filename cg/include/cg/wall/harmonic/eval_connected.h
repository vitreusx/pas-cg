#pragma once
#include "data.h"
#include <cg/simul/runtime.h>

namespace cg::wall::harmonic {
class eval_connected : public simul::sliceable_task {
public:
  real HH1;

public:
  vect::const_view<wall> walls;
  vect::const_view<connection> conns;
  vect::const_view<vec3r> r;
  vect::view<vec3r> F, wall_F;
  real *V;

public:
  void operator()() const;
  void omp_async() const;
  template <typename E> void iter(connection_expr<E> const &conn) const;

  void for_slice(int from, int to) const override;
  int total_size() const override;
  int slice_size() const override;
};
} // namespace cg::wall::harmonic