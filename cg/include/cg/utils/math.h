#pragma once
#include <cg/vect/vect.h>
#include <cmath>
#include <cuda.h>
#include <vcl/vectormath_exp.h>
#include <vcl/vectormath_trig.h>

#define C216 1.122462048309373
#define C216_INV 0.89089871814033927

namespace cg {
template <typename T>
__host__ __device__ auto isnan(T const &x) {
#ifndef __CUDA_ARCH__
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::is_nan(x);
  else
    return std::isnan(x);
#else
  return ::isnan(x);
#endif
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
__host__ __device__ auto select(Tc const &cond, T const &if_true,
                                T const &if_false) {
#ifndef __CUDA_ARCH__
  return select_impl<Tc, T>::impl(cond, if_true, if_false);
#else
  return cond ? if_true : if_false;
#endif
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
__host__ __device__ auto min(T1 const &x1, T2 const &x2) {
#ifndef __CUDA_ARCH__
  if constexpr (nitro::def::is_vcl_lane_v<T1>)
    return ::min(x1, x2);
  else
    return std::min(x1, x2);
#else
  return (x1 < x2) ? x1 : x2;
#endif
}

template <typename T1, typename T2>
__host__ __device__ auto max(T1 const &x1, T2 const &x2) {
#ifndef __CUDA_ARCH__
  if constexpr (nitro::def::is_vcl_lane_v<T1>)
    return ::max(x1, x2);
  else
    return std::max(x1, x2);
#else
  return (x1 > x2) ? x1 : x2;
#endif
}

template <typename Tv, typename Tmin, typename Tmax>
__host__ __device__ auto clamp(Tv value, Tmin min_value, Tmax max_value) {
  return cg::max(cg::min(value, max_value), min_value);
}

template <typename T>
__host__ __device__ auto sqrt(T const &x) {
#ifndef __CUDA_ARCH__
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::sqrt(x);
  else
    return std::sqrt(x);
#else
  return ::sqrt(x);
#endif
}

template <typename T>
__host__ __device__ auto rsqrt(T const &x) {
#ifndef __CUDA_ARCH__
  return (T)1 / cg::sqrt(x);
#else
  if constexpr (std::is_same_v<T, float>)
    return ::rsqrtf(x);
  else
    return (T)1 / cg::sqrt(x);
#endif
}

template <typename T>
__host__ __device__ auto apx_rsqrt(T const &x) {
#ifndef __CUDA_ARCH__
  if constexpr (nitro::def::is_vcl_lane_v<T> &&
                std::is_same_v<nitro::def::vcl_lane_type<T>, float>)
    return ::approx_rsqrt(x);
  else
    return cg::rsqrt(x);
#else
  return cg::rsqrt(x);
#endif
}

template <typename T>
__host__ __device__ auto apx_sqrt(T const &x) {
  return x * cg::apx_rsqrt(x);
}

template <typename T>
__host__ __device__ auto cbrt(T const &x) {
#ifndef __CUDA_ARCH__
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::cbrt(x);
  else
    return std::cbrt(x);
#else
  return ::cbrt(x);
#endif
}

template <size_t N, typename T>
__host__ __device__ auto ipow(T const &x) {
  if constexpr (N == 0)
    return (T)1;
  else
    return x * ipow<N - 1>(x);
}

template <typename T>
__host__ __device__ auto round(T x) {
#ifndef __CUDA_ARCH__
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::round(x);
  else
    return std::round(x);
#else
  return ::round(x);
#endif
}

template <typename T>
__host__ __device__ auto nearbyint(T x) {
#ifndef __CUDA_ARCH__
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return cg::round(x);
  else
    return std::nearbyint(x);
#else
  return cg::round(x);
#endif
}

template <typename T>
__host__ __device__ auto abs(T const &x) {
#ifndef __CUDA_ARCH__
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::abs(x);
  else
    return std::abs(x);
#else
  return ::abs(x);
#endif
}

template <typename T>
__host__ __device__ auto cos(T const &x) {
#ifndef __CUDA_ARCH__
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::cos(x);
  else
    return std::cos(x);
#else
  return ::abs(x);
#endif
}

template <typename T>
__host__ __device__ auto sin(T const &x) {
#ifndef __CUDA_ARCH__
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::sin(x);
  else
    return std::sin(x);
#else
  return ::sin(x);
#endif
}

template <typename T>
__host__ __device__ auto acos(T const &x) {
#ifndef __CUDA_ARCH__
  auto x_ = clamp(x, T(-1.0), T(1.0));
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::acos(x_);
  else
    return std::acos(x_);
#else
  return ::cos(x);
#endif
}

template <typename T>
__host__ __device__ auto exp(T const &x) {
#ifndef __CUDA_ARCH__
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::exp(x);
  else
    return std::exp(x);
#else
  return ::exp(x);
#endif
}

template <typename T>
__host__ __device__ auto ceil(T const &x) {
#ifndef __CUDA_ARCH__
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::ceil(x);
  else
    return std::ceil(x);
#else
  return ::ceil(x);
#endif
}

template <typename T>
__host__ __device__ auto floor(T const &x) {
#ifndef __CUDA_ARCH__
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::floor(x);
  else
    return std::floor(x);
#else
  return ::floor(x);
#endif
}

template <typename T>
__host__ __device__ auto log(T const &x) {
#ifndef __CUDA_ARCH__
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::log(x);
  else
    return std::log(x);
#else
  return ::log(x);
#endif
}

template <typename T>
__host__ __device__ auto isfinite(T const &x) {
#ifndef __CUDA_ARCH__
  if constexpr (nitro::def::is_vcl_lane_v<T>)
    return ::is_finite(x);
  else
    return std::isfinite(x);
#else
  return ::isfinite(x);
#endif
}

template <typename T>
__host__ __device__ auto sign(T const &x) {
#ifndef __CUDA_ARCH__
  if constexpr (nitro::def::is_vcl_lane_v<T>) {
    using IntLane = vect::lane<int, nitro::def::lane_size_v<T>,
                               nitro::def::lane_width_v<T>>;
    return IntLane(::sign_bit(x));
  } else
    return (T(0) < x) - (x < T(0));
#else
  return ::signbit(x);
#endif
}

} // namespace cg