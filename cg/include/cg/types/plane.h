#pragma once
#include "vec3.h"

namespace cg {
template <typename Scalar> class plane;

template <typename E> struct plane_expr {
  EXPR(origin, normal)

  template <typename F> auto projection(vec3_expr<F> const &v) const {
    return v - signed_dist(v) * normal();
  }

  template <typename F> auto signed_dist(vec3_expr<F> const &v) const {
    return dot(normal(), v - origin());
  }

  template <typename F> auto dist(vec3_expr<F> const &v) const {
    return abs(signed_dist(v));
  }
};

template <typename Scalar> class plane : public plane_expr<plane<Scalar>> {
public:
  INST(plane, FIELD(vec3<Scalar>, origin), FIELD(vec3<Scalar>, normal))
};
} // namespace cg