#pragma once
#include <nitro/nitro.h>

template <typename Scalar> class vec3;

template <typename E> struct vec3_expr : public ind_expr<E> {
  decltype(auto) x() { return static_cast<E &>(*this).x(); }
  decltype(auto) x() const { return static_cast<E const &>(*this).x(); }
  decltype(auto) y() { return static_cast<E &>(*this).y(); }
  decltype(auto) y() const { return static_cast<E const &>(*this).y(); }
  decltype(auto) z() { return static_cast<E &>(*this).z(); }
  decltype(auto) z() const { return static_cast<E const &>(*this).z(); }

  template <size_t I> decltype(auto) get() {
    if constexpr (I == 0)
      return x();
    else if constexpr (I == 1)
      return y();
    else if constexpr (I == 2)
      return z();
  }

  template <size_t I> decltype(auto) get() const {
    if constexpr (I == 0)
      return x();
    else if constexpr (I == 1)
      return y();
    else if constexpr (I == 2)
      return z();
  }

  template<typename F>
  vec3_expr<E>& omp_atomic_add(vec3_expr<F> const& f) {
#pragma omp atomic
    x() += f.x();
#pragma omp atomic
    y() += f.y();
#pragma omp atomic
    z() += f.z();

    return *this;
  }

  template<typename F>
  vec3_expr<E>& operator+=(vec3_expr<F> const& f) {
    x() += f.x();
    y() += f.y();
    z() += f.z();
    return *this;
  }

  template<typename F>
  vec3_expr<E>& operator-=(vec3_expr<F> const& f) {
    x() -= f.x();
    y() -= f.y();
    z() -= f.z();
    return *this;
  }

  template<typename F>
  vec3_expr<E>& operator*=(F const& f) {
    x() *= f;
    y() *= f;
    z() *= f;
    return *this;
  }

  template<typename F>
  vec3_expr<E>& operator/=(F const& f) {
    x() /= f;
    y() /= f;
    z() /= f;
    return this;
  }
};

template <typename E, typename Scalar> struct expr_impl<E, vec3<Scalar>> {
  using type = vec3_expr<E>;
};

template <typename E> struct vec3_auto_expr : public vec3_expr<E> {
  decltype(auto) x() { return static_cast<E &>(*this).template get<0>(); }
  decltype(auto) x() const {
    return static_cast<E const &>(*this).template get<0>();
  }
  decltype(auto) y() { return static_cast<E &>(*this).template get<1>(); }
  decltype(auto) y() const {
    return static_cast<E const &>(*this).template get<1>();
  }
  decltype(auto) z() { return static_cast<E &>(*this).template get<2>(); }
  decltype(auto) z() const {
    return static_cast<E const &>(*this).template get<2>();
  }

  template <size_t I> decltype(auto) get() {
    return static_cast<E &>(*this).template get<I>();
  }

  template <size_t I> decltype(auto) get() const {
    return static_cast<E const &>(*this).template get<I>();
  }
};

template <typename E, typename Scalar> struct auto_expr_impl<E, vec3<Scalar>> {
  using type = vec3_auto_expr<E>;
};

template <typename Scalar>
class vec3 : public vec3_auto_expr<vec3<Scalar>>,
             public tuple_wrapper<Scalar, Scalar, Scalar> {
public:
  using Base = tuple_wrapper<Scalar, Scalar, Scalar>;
  using Base::get;

  vec3() : Base((Scalar)0, (Scalar)0, (Scalar)0) {}

  vec3(Scalar &&x, Scalar &&y, Scalar &&z)
      : Base(std::forward<Scalar>(x), std::forward<Scalar>(y),
             std::forward<Scalar>(z)) {}

  template <typename E>
  vec3(vec3_expr<E> const &e) : vec3(e.x(), e.y(), e.z()) {}
};

template <typename Scalar>
struct is_indexed_impl<vec3<Scalar>> : std::true_type {};

using vec3f = vec3<float>;
using vec3d = vec3<float>;

#include "vec3_ops.h"