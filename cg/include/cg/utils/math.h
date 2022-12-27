#pragma once
#include <cg/vect/vect.h>
#include <cmath>
#include <vcl/vectormath_exp.h>
#include <vcl/vectormath_trig.h>

#define C216 1.122462048309373
#define C216_INV 0.89089871814033927

namespace cg {
template <typename T>
auto isnan(T const &x) {
  return std::isnan(x);
}

template <typename Tc, typename T>
auto select(Tc const &cond, T const &if_true, T const &if_false) {
  return cond ? if_true : if_false;
}

template <typename T>
auto clamp(T value, T min_value, T max_value) {
  return select(value < min_value, min_value,
                select(value < max_value, value, max_value));
}

template <typename T>
auto sqrt(T const &x) {
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::sqrt(x);
  else
    return std::sqrt(x);
}

template <typename T>
auto rsqrt(T const &x) {
  return (T)1 / sqrt(x);
}

template <typename T>
auto apx_rsqrt(T const &x) {
  return (T)1 / sqrt(x);
}

template <typename T>
auto apx_sqrt(T const &x) {
  return x * apx_rsqrt(x);
}

template <typename T>
auto cbrt(T const &x) {
  return std::cbrt(x);
}

template <size_t N, typename T>
auto ipow(T const &x) {
  if constexpr (N == 0)
    return (T)1;
  else
    return x * ipow<N - 1>(x);
}

template <typename T>
auto nearbyint(T x) {
  return std::nearbyint(x);
}

template <typename T>
auto abs(T const &x) {
  return std::abs(x);
}

template <typename T>
auto cos(T const &x) {
  return std::cos(x);
}

template <typename T>
auto sin(T const &x) {
  return std::sin(x);
}

template <typename T>
auto acos(T const &x) {
  return std::acos(clamp(x, T(-1.0), T(1.0)));
}

template <typename T>
auto exp(T const &x) {
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::exp(x);
  else
    return std::exp(x);
}

template <typename T1, typename T2>
auto min(T1 const &x1, T2 const &x2) {
  return std::min(x1, x2);
}

template <typename T1, typename T2>
auto max(T1 const &x1, T2 const &x2) {
  return std::max(x1, x2);
}

template <typename T>
auto ceil(T const &x) {
  return std::ceil(x);
}

template <typename T>
auto floor(T const &x) {
  return std::floor(x);
}

template <typename T>
auto log(T const &x) {
  return std::log(x);
}

template <typename T>
auto isfinite(T const &x) {
  return std::isfinite(x);
}

template <typename T>
int sign(T const &value) {
  return (T(0) < value) - (value < T(0));
}

} // namespace cg