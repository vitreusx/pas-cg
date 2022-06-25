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
auto operator+(vec3_expr<E1> const &e1, vec3_expr<E2> const &e2) {
  auto x = e1.x() + e2.x();
  auto y = e1.y() + e2.y();
  auto z = e1.z() + e2.z();
  return vec3<decltype(x)>(x, y, z);
}

template <typename E1, typename E2>
auto operator-(vec3_expr<E1> const &e1, vec3_expr<E2> const &e2) {
  auto x = e1.x() - e2.x();
  auto y = e1.y() - e2.y();
  auto z = e1.z() - e2.z();
  return vec3<decltype(x)>(x, y, z);
}

template <typename E> auto operator-(vec3_expr<E> const &e) {
  auto x = -e.x();
  auto y = -e.y();
  auto z = -e.z();
  return vec3<decltype(x)>(x, y, z);
}

template <typename S, typename E>
auto operator*(S const &s, vec3_expr<E> const &e) {
  auto x = s * e.x();
  auto y = s * e.y();
  auto z = s * e.z();
  return vec3<decltype(x)>(x, y, z);
}

template <typename E, typename S>
auto operator*(vec3_expr<E> const &e, S const &s) {
  auto x = e.x() * s;
  auto y = e.y() * s;
  auto z = e.z() * s;
  return vec3<decltype(x)>(x, y, z);
}

template <typename E, typename S>
auto operator/(vec3_expr<E> const &e, S const &s) {
  auto x = e.x() / s;
  auto y = e.y() / s;
  auto z = e.z() / s;
  return vec3<decltype(x)>(x, y, z);
}

template <typename E1, typename E2>
auto cross(vec3_expr<E1> const &e1, vec3_expr<E2> const &e2) {
  auto x = e1.y() * e2.z() - e1.z() * e2.y();
  auto y = e1.z() * e2.x() - e1.x() * e2.z();
  auto z = e1.x() * e2.y() - e1.y() * e2.x();
  return vec3<decltype(x)>(x, y, z);
}

template <typename E> auto unit(vec3_expr<E> const &e) {
  return norm_inv(e) * e;
}

template <typename S, typename E>
Eigen::Vector3<S> convert(vec3_expr<E> const &e) {
  return {e.x(), e.y(), e.z()};
}
} // namespace cg