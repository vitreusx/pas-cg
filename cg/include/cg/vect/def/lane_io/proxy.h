#pragma once
#include "../lane/proxy.h"
#include "vcl.h"

namespace nitro::def {
template <typename Lane, typename T, std::size_t N, std::size_t W,
          typename =
              std::enable_if_t<std::is_same_v<Lane, vcl_lane<repr_t<T>, N, W>>>>
void load(Lane &data, T const *src) {
  load(data, reinterpret_cast<repr_t<T> const *>(src));
}

template <typename Lane, typename T, std::size_t N, std::size_t W,
          typename =
              std::enable_if_t<std::is_same_v<Lane, vcl_lane<repr_t<T>, N, W>>>>
void store(Lane const &data, T *dst) {
  store(data, reinterpret_cast<repr_t<T> *>(dst));
}

template <typename Lane, typename T, std::size_t N, std::size_t W,
          typename =
              std::enable_if_t<std::is_same_v<Lane, vcl_lane<repr_t<T>, N, W>>>>
Lane construct(T const *src) {
  return construct<Lane>(reinterpret_cast<repr_t<T> const *>(src));
}
} // namespace nitro::def