#pragma once
#include "lane/array.h"
#include "lane/proxy.h"
#include "lane/vcl.h"
#include "lane_io/array.h"
#include "lane_io/proxy.h"
#include "lane_io/vcl.h"

namespace nitro::def {
template <typename T, std::size_t N, std::size_t W>
struct _lane {
  template <typename U1, typename = void>
  struct vcl_or {
    template <typename U2, typename = void>
    struct proxy_or {
      using type = array_lane<U2, N>;
    };

    template <typename U2>
    struct proxy_or<U2,
                    std::void_t<typename _vcl_lane<repr_t<U2>, N, W>::type>> {
      using type = vcl_lane<repr_t<U2>, N, W>;
    };

    using type = typename proxy_or<U1>::type;
  };

  template <typename U1>
  struct vcl_or<U1, std::void_t<typename _vcl_lane<U1, N, W>::type>> {
    using type = vcl_lane<U1, N, W>;
  };

  using type = typename vcl_or<T>::type;
};

#ifdef __AVX512F__
inline constexpr std::size_t opt_width_v = 512;
#else
inline constexpr std::size_t opt_width_v = 256;
#endif

template <typename T, std::size_t N, std::size_t W = opt_width_v>
using lane = typename _lane<T, N, W>::type;

template <typename Lane, typename T>
void load(Lane &data, T const *src) {
  if constexpr (has_repr_v<T>)
    proxy_load(data, src);
  else if constexpr (is_vcl_lane_v<Lane>)
    vcl_load(data, src);
  else
    array_load(data, src);
}

template <typename Lane, typename T>
void store(Lane const &data, T *dst) {
  if constexpr (has_repr_v<T>)
    proxy_store(data, dst);
  else if constexpr (is_vcl_lane_v<Lane>)
    vcl_store(data, dst);
  else
    array_store(data, dst);
}

template <typename Lane, typename T>
Lane construct(T const *src) {
  if constexpr (has_repr_v<T>)
    return proxy_construct<Lane>(src);
  else if constexpr (is_vcl_lane_v<Lane>)
    return vcl_construct<Lane>(src);
  else
    return array_construct<Lane>(src);
}

template <typename Lane, typename T, typename Idxes>
Lane gather(T const *src, Idxes const &idxes) {
  if constexpr (has_repr_v<T>)
    return proxy_gather<Lane>(src, idxes);
  else if constexpr (is_vcl_lane_v<Lane>)
    return vcl_gather<Lane>(src, idxes);
  else
    return array_gather<Lane>(src, idxes);
}

template <typename Lane, typename T, typename Idxes>
void scatter(Lane const &data, T *dst, Idxes const &idxes) {
  if constexpr (has_repr_v<T>)
    proxy_scatter(data, dst, idxes);
  else if constexpr (is_vcl_lane_v<Lane>)
    vcl_scatter(data, dst, idxes);
  else
    array_scatter(data, dst, idxes);
}

template <typename Lane, typename T, typename Idxes, typename Mask>
Lane masked_gather(T const *src, Idxes const &idxes, Mask const &mask) {
  if constexpr (has_repr_v<T>)
    return proxy_masked_gather<Lane>(src, idxes, mask);
  else if constexpr (is_vcl_lane_v<Lane>)
    return vcl_masked_gather<Lane>(src, idxes, mask);
  //  else
  //    return array_masked_gather<Lane>(src, idxes, mask);
}

template <typename Lane, typename T, typename Idxes, typename Mask>
void masked_scatter(Lane const &data, T *dst, Idxes const &idxes,
                    Mask const &mask) {
  if constexpr (has_repr_v<T>)
    proxy_masked_scatter(data, dst, idxes, mask);
  else if constexpr (is_vcl_lane_v<Lane>)
    vcl_masked_scatter(data, dst, idxes, mask);
  //  else
  //    array_masked_scatter(dst, idxes, mask);
}

template<typename Lane, typename T, typename Idx>
Lane lookup(T const *src, Idx const *idxes) {
  if constexpr (has_repr_v<T>)
    return proxy_lookup<Lane>(src, idxes);
  else if constexpr (is_vcl_lane_v<Lane>)
    return vcl_lookup<Lane>(src, idxes);
  else
    return array_lookup<Lane>(src, idxes);
}

} // namespace nitro::def