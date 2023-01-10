#pragma once
#include "vec3_def.h"
#include <cg/utils/math.h>

namespace cg {
template <typename E>
__host__ __device__ auto norm(vec3_expr<E> const &e) {
  return cg::sqrt(norm_squared(e));
}

template <typename E>
__host__ __device__ auto norm_squared(vec3_expr<E> const &e) {
  return e.x() * e.x() + e.y() * e.y() + e.z() * e.z();
}

template <typename E>
__host__ __device__ auto l1_norm(vec3_expr<E> const &e) {
  return cg::abs(e.x()) + cg::abs(e.y()) + cg::abs(e.z());
}

template <typename E1, typename E2>
__host__ __device__ auto dot(vec3_expr<E1> const &e1, vec3_expr<E2> const &e2) {
  return e1.x() * e2.x() + e1.y() * e2.y() + e1.z() * e2.z();
}

template <typename E>
__host__ __device__ auto norm_inv(vec3_expr<E> const &e) {
  return cg::rsqrt(norm_squared(e));
}

template <typename E>
__host__ __device__ auto apx_norm(vec3_expr<E> const &e) {
  return cg::apx_sqrt(norm_squared(e));
}

template <typename E>
__host__ __device__ auto apx_norm_inv(vec3_expr<E> const &e) {
  return cg::apx_rsqrt(norm_squared(e));
}

template <typename E1, typename E2>
class sum_expr : public vec3_expr<sum_expr<E1, E2>> {
public:
  __host__ __device__ sum_expr(E1 const &e1, E2 const &e2) : e1{e1}, e2{e2} {};

  __host__ __device__ auto x() const {
    return e1.x() + e2.x();
  }

  __host__ __device__ auto y() const {
    return e1.y() + e2.y();
  }

  __host__ __device__ auto z() const {
    return e1.z() + e2.z();
  }

private:
  E1 e1;
  E2 e2;
};

template <typename E1, typename E2>
__host__ __device__ auto operator+(vec3_expr<E1> const &e1,
                                   vec3_expr<E2> const &e2) {
  return sum_expr<E1, E2>(static_cast<E1 const &>(e1),
                          static_cast<E2 const &>(e2));
}

template <typename E1, typename E2>
class diff_expr : public vec3_expr<diff_expr<E1, E2>> {
public:
  __host__ __device__ diff_expr(E1 const &e1, E2 const &e2) : e1{e1}, e2{e2} {};

  __host__ __device__ auto x() const {
    return e1.x() - e2.x();
  }

  __host__ __device__ auto y() const {
    return e1.y() - e2.y();
  }

  __host__ __device__ auto z() const {
    return e1.z() - e2.z();
  }

private:
  E1 e1;
  E2 e2;
};

template <typename E1, typename E2>
__host__ __device__ auto operator-(vec3_expr<E1> const &e1,
                                   vec3_expr<E2> const &e2) {
  return diff_expr<E1, E2>(static_cast<E1 const &>(e1),
                           static_cast<E2 const &>(e2));
}

template <typename E>
class neg_expr : public vec3_expr<neg_expr<E>> {
public:
  __host__ __device__ explicit neg_expr(E const &e) : e{e} {};

  __host__ __device__ auto x() const {
    return -e.x();
  }

  __host__ __device__ auto y() const {
    return -e.y();
  }

  __host__ __device__ auto z() const {
    return -e.z();
  }

private:
  E e;
};

template <typename E>
__host__ __device__ auto operator-(vec3_expr<E> const &e) {
  return neg_expr<E>(static_cast<E const &>(e));
}

template <typename S, typename E>
class scalar_lmul_expr : public vec3_expr<scalar_lmul_expr<S, E>> {
public:
  __host__ __device__ explicit scalar_lmul_expr(S const &s, E const &e)
      : s{s}, e{e} {};

  __host__ __device__ auto x() const {
    return s * e.x();
  }

  __host__ __device__ auto y() const {
    return s * e.y();
  }

  __host__ __device__ auto z() const {
    return s * e.z();
  }

private:
  S s;
  E e;
};

template <typename S, typename E>
__host__ __device__ auto operator*(S const &s, vec3_expr<E> const &e) {
  return scalar_lmul_expr<S, E>(s, static_cast<E const &>(e));
}

template <typename E, typename S>
class scalar_rmul_expr : public vec3_expr<scalar_rmul_expr<E, S>> {
public:
  __host__ __device__ explicit scalar_rmul_expr(E const &e, S const &s)
      : e{e}, s{s} {};

  __host__ __device__ auto x() const {
    return e.x() * s;
  }

  __host__ __device__ auto y() const {
    return e.y() * s;
  }

  __host__ __device__ auto z() const {
    return e.z() * s;
  }

private:
  E e;
  S s;
};

template <typename E, typename S>
__host__ __device__ auto operator*(vec3_expr<E> const &e, S const &s) {
  return scalar_rmul_expr<E, S>(static_cast<E const &>(e), s);
}

template <typename E, typename S>
class scalar_div_expr : public vec3_expr<scalar_div_expr<E, S>> {
public:
  __host__ __device__ explicit scalar_div_expr(E const &e, S const &s)
      : e{e}, s{s} {};

  __host__ __device__ auto x() const {
    return e.x() / s;
  }

  __host__ __device__ auto y() const {
    return e.y() / s;
  }

  __host__ __device__ auto z() const {
    return e.z() / s;
  }

private:
  E e;
  S s;
};

template <typename E, typename S>
__host__ __device__ auto operator/(vec3_expr<E> const &e, S const &s) {
  return scalar_div_expr<E, S>(static_cast<E const &>(e), s);
}

template <typename E1, typename E2>
class cross_expr : public vec3_expr<cross_expr<E1, E2>> {
public:
  __host__ __device__ explicit cross_expr(E1 const &e1, E2 const &e2)
      : e1{e1}, e2{e2} {};

  __host__ __device__ auto x() const {
    return e1.y() * e2.z() - e1.z() * e2.y();
  }

  __host__ __device__ auto y() const {
    return e1.z() * e2.x() - e1.x() * e2.z();
  }

  __host__ __device__ auto z() const {
    return e1.x() * e2.y() - e1.y() * e2.x();
  }

private:
  E1 e1;
  E2 e2;
};

template <typename E1, typename E2>
__host__ __device__ auto cross(vec3_expr<E1> const &e1,
                               vec3_expr<E2> const &e2) {
  return cross_expr<E1, E2>(static_cast<E1 const &>(e1),
                            static_cast<E2 const &>(e2));
}

template <typename E>
__host__ __device__ auto unit(vec3_expr<E> const &e) {
  return norm_inv(e) * e;
}

template <typename E1, typename E2>
__host__ __device__ auto cast(vec3_expr<E1> const &v,
                              vec3_expr<E2> const &onto) {
  return onto * dot(v, onto);
}

template <typename S, typename E>
Eigen::Vector3<S> convert(vec3_expr<E> const &e) {
  return {e.x(), e.y(), e.z()};
}

} // namespace cg
