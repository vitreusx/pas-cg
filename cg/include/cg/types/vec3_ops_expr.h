#pragma once
#include "vec3_def.h"
#include <cg/utils/math.h>

namespace cg {
template <typename E> auto norm(vec3_expr<E> const &e) {
  return cg::sqrt(norm_squared(e));
}

template <typename E> auto norm_squared(vec3_expr<E> const &e) {
  return e.x() * e.x() + e.y() * e.y() + e.z() * e.z();
}

template <typename E> auto l1_norm(vec3_expr<E> const &e) {
  return cg::abs(e.x()) + cg::abs(e.y()) + cg::abs(e.z());
}

template <typename E1, typename E2>
auto dot(vec3_expr<E1> const &e1, vec3_expr<E2> const &e2) {
  return e1.x() * e2.x() + e1.y() * e2.y() + e1.z() * e2.z();
}

template <typename E> auto norm_inv(vec3_expr<E> const &e) {
  return cg::rsqrt(norm_squared(e));
}

template <typename E> auto apx_norm(vec3_expr<E> const &e) {
  return cg::apx_sqrt(norm_squared(e));
}

template <typename E> auto apx_norm_inv(vec3_expr<E> const &e) {
  return cg::apx_rsqrt(norm_squared(e));
}

template <typename E1, typename E2>
class sum_expr : public vec3_expr<sum_expr<E1, E2>> {
public:
  sum_expr(E1 const &e1, E2 const &e2) : e1{e1}, e2{e2} {};

  auto x() const {
    return e1.x() + e2.x();
  }

  auto y() const {
    return e1.y() + e2.y();
  }

  auto z() const {
    return e1.z() + e2.z();
  }

private:
  E1 e1;
  E2 e2;
};

template <typename E1, typename E2>
auto operator+(vec3_expr<E1> const &e1, vec3_expr<E2> const &e2) {
  return sum_expr<E1, E2>(static_cast<E1 const &>(e1),
                          static_cast<E2 const &>(e2));
}

template <typename E1, typename E2>
class diff_expr : public vec3_expr<diff_expr<E1, E2>> {
public:
  diff_expr(E1 const &e1, E2 const &e2) : e1{e1}, e2{e2} {};

  auto x() const {
    return e1.x() - e2.x();
  }

  auto y() const {
    return e1.y() - e2.y();
  }

  auto z() const {
    return e1.z() - e2.z();
  }

private:
  E1 e1;
  E2 e2;
};

template <typename E1, typename E2>
auto operator-(vec3_expr<E1> const &e1, vec3_expr<E2> const &e2) {
  return diff_expr<E1, E2>(static_cast<E1 const &>(e1),
                           static_cast<E2 const &>(e2));
}

template <typename E> class neg_expr : public vec3_expr<neg_expr<E>> {
public:
  explicit neg_expr(E const &e) : e{e} {};

  auto x() const {
    return -e.x();
  }

  auto y() const {
    return -e.y();
  }

  auto z() const {
    return -e.z();
  }

private:
  E e;
};

template <typename E> auto operator-(vec3_expr<E> const &e) {
  return neg_expr<E>(static_cast<E const &>(e));
}

template <typename S, typename E>
class scalar_lmul_expr : public vec3_expr<scalar_lmul_expr<S, E>> {
public:
  explicit scalar_lmul_expr(S const &s, E const &e) : s{s}, e{e} {};

  auto x() const {
    return s * e.x();
  }

  auto y() const {
    return s * e.y();
  }

  auto z() const {
    return s * e.z();
  }

private:
  S s;
  E e;
};

template <typename S, typename E>
auto operator*(S const &s, vec3_expr<E> const &e) {
  return scalar_lmul_expr<S, E>(s, static_cast<E const &>(e));
}

template <typename E, typename S>
class scalar_rmul_expr : public vec3_expr<scalar_rmul_expr<E, S>> {
public:
  explicit scalar_rmul_expr(E const &e, S const &s) : e{e}, s{s} {};

  auto x() const {
    return e.x() * s;
  }

  auto y() const {
    return e.y() * s;
  }

  auto z() const {
    return e.z() * s;
  }

private:
  E e;
  S s;
};

template <typename E, typename S>
auto operator*(vec3_expr<E> const &e, S const &s) {
  return scalar_rmul_expr<E, S>(static_cast<E const &>(e), s);
}

template <typename E, typename S>
class scalar_div_expr : public vec3_expr<scalar_div_expr<E, S>> {
public:
  explicit scalar_div_expr(E const &e, S const &s) : e{e}, s{s} {};

  auto x() const {
    return e.x() / s;
  }

  auto y() const {
    return e.y() / s;
  }

  auto z() const {
    return e.z() / s;
  }

private:
  E e;
  S s;
};

template <typename E, typename S>
auto operator/(vec3_expr<E> const &e, S const &s) {
  return scalar_div_expr<E, S>(static_cast<E const &>(e), s);
}

template <typename E1, typename E2>
class cross_expr : public vec3_expr<cross_expr<E1, E2>> {
public:
  explicit cross_expr(E1 const &e1, E2 const &e2) : e1{e1}, e2{e2} {};

  auto x() const {
    return e1.y() * e2.z() - e1.z() * e2.y();
  }

  auto y() const {
    return e1.z() * e2.x() - e1.x() * e2.z();
  }

  auto z() const {
    return e1.x() * e2.y() - e1.y() * e2.x();
  }

private:
  E1 e1;
  E2 e2;
};

template <typename E1, typename E2>
auto cross(vec3_expr<E1> const &e1, vec3_expr<E2> const &e2) {
  return cross_expr<E1, E2>(static_cast<E1 const &>(e1),
                            static_cast<E2 const &>(e2));
}

template <typename E> auto unit(vec3_expr<E> const &e) {
  return norm_inv(e) * e;
}

template <typename E1, typename E2>
auto cast(vec3_expr<E1> const &v, vec3_expr<E2> const &onto) {
  return onto * dot(v, onto);
}

template <typename S, typename E>
Eigen::Vector3<S> convert(vec3_expr<E> const &e) {
  return {e.x(), e.y(), e.z()};
}

} // namespace cg
