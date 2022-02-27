#pragma once
#include <cmath>
#include <nitro/nitro.h>

namespace cg {
template <typename T> auto sqrt(T const &x) { return std::sqrt(x); }

template <typename T, size_t N> auto sqrt(nitro::lane<T, N> const &x) {
  return ::square(x);
}

template <typename T> auto rsqrt(T const &x) { return T(1) / sqrt(x); }

template <typename T> auto apx_rsqrt(T const &x) { return (T)1 / sqrt(x); }

template <typename T, size_t N> auto apx_rsqrt(nitro::lane<T, N> const &x) {
  return apx_rsqrt(x);
}

template <typename T> auto apx_sqrt(T const &x) { return x * apx_rsqrt(x); }

template <typename T> auto cbrt(T const &x) { return std::cbrt(x); }

template <typename T, size_t N> auto cbrt(nitro::lane<T, N> const &x) {
  return ::cbrt(x);
}

template <size_t N, typename T> auto ipow(T const &x) {
  if constexpr (N == 0)
    return (T)1;
  else
    return x * ipow<N - 1>(x);
}

template <size_t N, typename T, size_t M>
auto ipow(nitro::lane<T, M> const &x) {
  return ::pow_const(x, N);
}

template <typename T> auto nearbyint(T x) { return std::nearbyint(x); }

template <typename T, size_t N> auto nearbyint(nitro::lane<T, N> const &x) {
  return ::round(x);
}

template <typename T> auto abs(T const &x) { return std::abs(x); }

template <typename T, size_t N> auto abs(nitro::lane<T, N> const &x) {
  return ::abs(x);
}

template <typename T> auto cos(T const &x) { return std::cos(x); }

template <typename T, size_t N> auto cos(nitro::lane<T, N> const &x) {
  return ::cos(x);
}

template <typename T> auto sin(T const &x) { return std::sin(x); }

template <typename T, size_t N> auto sin(nitro::lane<T, N> const &x) {
  return ::sin(x);
}

template <typename T> auto acos(T const &x) {
  auto x_ = x < T(-1.0) ? T(-1.0) : (x > T(1.0) ? T(1.0) : x);
  return std::acos(x_);
}

template <typename T, size_t N> auto acos(nitro::lane<T, N> const &x) {
  return ::acos(x);
}

template <typename T> auto exp(T const &x) { return std::exp(x); }

template <typename T, size_t N> auto exp(nitro::lane<T, N> const &x) {
  return ::exp(x);
}

template <typename T1, typename T2> auto min(T1 const &x1, T2 const &x2) {
  return std::min(x1, x2);
}

template <typename T, size_t N>
auto min(nitro::lane<T, N> const &x1, nitro::lane<T, N> const &x2) {
  return ::min(x1, x2);
}

template <typename T1, typename T2> auto max(T1 const &x1, T2 const &x2) {
  return std::max(x1, x2);
}

template <typename T, size_t N>
auto max(nitro::lane<T, N> const &x1, nitro::lane<T, N> const &x2) {
  return ::max(x1, x2);
}

template <typename T> auto ceil(T const &x) { return std::ceil(x); }

template <typename T, size_t N> auto ceil(nitro::lane<T, N> const &x) {
  return ::ceil(x);
}

template <typename T> auto floor(T const &x) { return std::floor(x); }

template <typename T, size_t N> auto floor(nitro::lane<T, N> const &x) {
  return ::floor(x);
}

template <typename T> auto log(T const &x) { return std::log(x); }

template <typename T, size_t N> auto log(nitro::lane<T, N> const &x) {
  return ::log(x);
}

template <typename T> auto isnan(T const &x) { return std::isnan(x); }

template <typename T, size_t N> auto isnan(nitro::lane<T, N> const &x) {
  return ::is_nan(x);
}

template <typename Tc, typename T>
auto select(Tc const &cond, T const &if_true, T const &if_false) {
  return cond ? if_true : if_false;
}

template <typename T, size_t N>
auto select(nitro::mask<nitro::lane<T, N>> const &cond,
            nitro::lane<T, N> const &if_true,
            nitro::lane<T, N> const &if_false) {
  return ::select(cond, if_true, if_false);
}

template <typename T> auto clamp(T value, T min_value, T max_value) {
  return select(value < min_value, min_value,
                select(value < max_value, value, max_value));
}
} // namespace cg