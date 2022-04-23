#pragma once
#include <Eigen/Core>
#include <cg/vect/vect.h>

namespace cg {
template <typename Scalar> class vec3;

template <typename E> struct vec3_expr : public nitro::ind_expr<E> {
  EXPR_BODY(x, y, z)

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
    return this;
  }
};

template <typename E> struct vec3_auto_expr : public vec3_expr<E> {
  AUTO_EXPR_BODY(x, y, z)
};

template <typename Scalar>
using vec3_base = nitro::tuple_wrapper<Scalar, Scalar, Scalar>;

template <typename Scalar>
class vec3 : public vec3_auto_expr<vec3<Scalar>>, public vec3_base<Scalar> {
public:
  using Base = vec3_base<Scalar>;
  using Base::get;

  static vec3<Scalar> Zero() { return {(Scalar)0, (Scalar)0, (Scalar)0}; }

  static vec3<Scalar> UnitX() { return {(Scalar)1, (Scalar)0, (Scalar)0}; }

  static vec3<Scalar> UnitY() { return {(Scalar)0, (Scalar)1, (Scalar)0}; }

  static vec3<Scalar> UnitZ() { return {(Scalar)0, (Scalar)0, (Scalar)1}; }

  vec3() : Base((Scalar)0, (Scalar)0, (Scalar)0) {}

  vec3(Scalar x, Scalar y, Scalar z)
      : Base(std::forward<Scalar>(x), std::forward<Scalar>(y),
             std::forward<Scalar>(z)) {}

  template <typename E>
  vec3(vec3_expr<E> const &e) : vec3(e.x(), e.y(), e.z()) {}

  template <typename T>
  vec3(Eigen::Vector3<T> const &eigen)
      : vec3((T)eigen.x(), (T)eigen.y(), (T)eigen.z()) {}
};

using vec3f = vec3<float>;
using vec3d = vec3<double>;

} // namespace cg

namespace nitro {

template <typename Scalar>
struct is_indexed_impl<cg::vec3<Scalar>> : std::true_type {};

template <typename E, typename Scalar> struct expr_impl<E, cg::vec3<Scalar>> {
  using type = cg::vec3_expr<E>;
};

template <typename E, typename Scalar>
struct auto_expr_impl<E, cg::vec3<Scalar>> {
  using type = cg::vec3_auto_expr<E>;
};

} // namespace nitro
