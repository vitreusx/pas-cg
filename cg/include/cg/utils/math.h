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
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::is_nan(x);
  else
    return std::isnan(x);
}

template <typename Tc, typename T, typename = void>
struct select_impl {
  static auto impl(Tc const &cond, T const &if_true, T const &if_false) {
    if constexpr (nitro::def::is_vcl_lane_v<Tc>) {
      using MaskT = decltype(std::declval<T>() < std::declval<T>());
      return ::select((MaskT)cond, if_true, if_false);
    } else {
      return cond ? if_true : if_false;
    }
  }
};

template <typename Tc, typename T>
auto select(Tc const &cond, T const &if_true, T const &if_false) {
  return select_impl<Tc, T>::impl(cond, if_true, if_false);
}

template <typename Tc, typename T>
struct select_impl<Tc, T, std::enable_if_t<vect::is_indexed_v<T>>> {
  template <std::size_t... Idxes>
  static auto aux(Tc const &cond, T const &if_true, T const &if_false,
                  vect::ind_seq<Idxes...>) {
    return T(cg::select(cond, if_true.template get<Idxes>(),
                        if_false.template get<Idxes>())...);
  }

  static auto impl(Tc const &cond, T const &if_true, T const &if_false) {
    return aux(cond, if_true, if_false, vect::idxes_t<T>{});
  }
};

template <typename T1, typename T2>
auto min(T1 const &x1, T2 const &x2) {
  if constexpr (nitro::def::is_vcl_lane_v<T1>)
    return ::min(x1, x2);
  else
    return std::min(x1, x2);
}

template <typename T1, typename T2>
auto max(T1 const &x1, T2 const &x2) {
  if constexpr (nitro::def::is_vcl_lane_v<T1>)
    return ::max(x1, x2);
  else
    return std::max(x1, x2);
}

template <typename Tv, typename Tmin, typename Tmax>
auto clamp(Tv value, Tmin min_value, Tmax max_value) {
  return cg::max(cg::min(value, max_value), min_value);
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
  return (T)1 / cg::sqrt(x);
}

template <typename T>
auto apx_rsqrt(T const &x) {
  if constexpr (nitro::def::is_vcl_lane_v<T> &&
                std::is_same_v<nitro::def::vcl_lane_type<T>, float>)
    return ::approx_rsqrt(x);
  else
    return (T)1 / sqrt(x);
}

template <typename T>
auto apx_sqrt(T const &x) {
  return x * cg::apx_rsqrt(x);
}

template <typename T>
auto cbrt(T const &x) {
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::cbrt(x);
  else
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
auto round(T x) {
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::round(x);
  else
    return std::round(x);
}

template <typename T>
auto nearbyint(T x) {
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return cg::round(x);
  else
    return std::nearbyint(x);
}

template <typename T>
auto abs(T const &x) {
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::abs(x);
  else
    return std::abs(x);
}

template <typename T>
auto cos(T const &x) {
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::cos(x);
  else
    return std::cos(x);
}

template <typename T>
auto sin(T const &x) {
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::sin(x);
  else
    return std::sin(x);
}

template <typename T>
auto acos(T const &x) {
  auto x_ = clamp(x, T(-1.0), T(1.0));
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::acos(x_);
  else
    return std::acos(x_);
}

template <typename T>
auto exp(T const &x) {
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::exp(x);
  else
    return std::exp(x);
}

template <typename T>
auto ceil(T const &x) {
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::ceil(x);
  else
    return std::ceil(x);
}

template <typename T>
auto floor(T const &x) {
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::floor(x);
  else
    return std::floor(x);
}

template <typename T>
auto log(T const &x) {
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::log(x);
  else
    return std::log(x);
}

template <typename T>
auto isfinite(T const &x) {
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::is_finite(x);
  else
    return std::isfinite(x);
}

template <typename T>
auto sign(T const &x) {
  if constexpr (nitro::def::is_vcl_lane_v<T>) {
    using IntLane = vect::lane<int, nitro::def::lane_size_v<T>,
                               nitro::def::lane_width_v<T>>;
    return IntLane(::sign_bit(x));
  } else
    return (T(0) < x) - (x < T(0));
}

} // namespace cg