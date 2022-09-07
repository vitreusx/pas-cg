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
} // namespace nitro::def