#pragma once
#include <Eigen/Core>
#include <cg/vect/vect.h>

namespace cg {
template <typename Scalar> class vec3;

template <typename E> struct vec3_expr {
  EXPR(x, y, z)

  template <typename F> vec3_expr<E> &omp_atomic_add(vec3_expr<F> const &f) {
#pragma omp atomic
    x() += f.x();
#pragma omp atomic
    y() += f.y();
#pragma omp atomic
    z() += f.z();

    return *this;
  }

  template <typename F> vec3_expr<E> &operator+=(vec3_expr<F> const &f) {
    x() += f.x();
    y() += f.y();
    z() += f.z();
    return *this;
  }

  template <typename F> vec3_expr<E> &operator-=(vec3_expr<F> const &f) {
    x() -= f.x();
    y() -= f.y();
    z() -= f.z();
    return *this;
  }

  template <typename F> vec3_expr<E> &operator*=(F const &f) {
    x() *= f;
    y() *= f;
    z() *= f;
    return *this;
  }

  template <typename F> vec3_expr<E> &operator/=(F const &f) {
    x() /= f;
    y() /= f;
    z() /= f;
    return *this;
  }
};

template <typename Scalar> class vec3 : public vec3_expr<vec3<Scalar>> {
public:
  INST(vec3, FIELD(Scalar, x), FIELD(Scalar, y), FIELD(Scalar, z))

  static vec3<Scalar> Zero() {
    return {(Scalar)0, (Scalar)0, (Scalar)0};
  }

  static vec3<Scalar> UnitX() {
    return {(Scalar)1, (Scalar)0, (Scalar)0};
  }

  static vec3<Scalar> UnitY() {
    return {(Scalar)0, (Scalar)1, (Scalar)0};
  }

  static vec3<Scalar> UnitZ() {
    return {(Scalar)0, (Scalar)0, (Scalar)1};
  }

  vec3() : vec3((Scalar)0.0, (Scalar)0.0, (Scalar)0.0) {}

  template <typename T>
  vec3(Eigen::Vector3<T> const &eigen)
      : vec3((T)eigen.x(), (T)eigen.y(), (T)eigen.z()) {}
};

using vec3f = vec3<float>;
using vec3d = vec3<double>;

} // namespace cg
