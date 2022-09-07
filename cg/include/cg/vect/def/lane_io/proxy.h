#pragma once
#include "../lane/proxy.h"
#include "vcl.h"

namespace nitro::def {
template <typename Lane, typename T>
void proxy_load(Lane &data, T const *src) {
  vcl_load(data, reinterpret_cast<repr_t<T> const *>(src));
}

template <typename Lane, typename T>
void proxy_store(Lane const &data, T *dst) {
  store(data, reinterpret_cast<repr_t<T> *>(dst));
}

template <typename Lane, typename T>
Lane proxy_construct(T const *src) {
  return vcl_construct<Lane>(reinterpret_cast<repr_t<T> const *>(src));
}

template <typename Lane, typename T, typename Idxes>
Lane proxy_gather(T const *src, Idxes const &idxes) {
  return vcl_gather<Lane>(reinterpret_cast<repr_t<T> const *>(src), idxes);
}

template <typename Lane, typename T, typename Idxes, typename Mask>
Lane proxy_masked_gather(T const *src, Idxes const &idxes, Mask const &mask) {
  return vcl_masked_gather<Lane>(reinterpret_cast<repr_t<T> const *>(src),
                                 idxes, mask);
}

template <typename Lane, typename T, typename Idxes>
void proxy_scatter(Lane const &data, T *dst, Idxes const &idxes) {
  vcl_scatter(data, reinterpret_cast<repr_t<T> *>(dst), idxes);
}

template <typename Lane, typename T, typename Idxes, typename Mask>
void proxy_masked_scatter(Lane const &data, T *dst, Idxes const &idxes,
                          Mask const &mask) {
  vcl_masked_scatter(data, reinterpret_cast<repr_t<T> *>(dst), idxes, mask);
}
} // namespace nitro::def