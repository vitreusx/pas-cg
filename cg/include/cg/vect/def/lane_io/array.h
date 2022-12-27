#pragma once
#include "../lane/array.h"
#include "decl.h"
#include "vcl.h"

namespace nitro::def {

template <typename T, std::size_t N>
void array_load(std::array<T, N> &data, T const *src) {
  for (int idx = 0; idx < N; ++idx)
    data[idx] = src[idx];
}

template <typename T, std::size_t N>
void array_store(std::array<T, N> const &data, T *dst) {
  for (int idx = 0; idx < N; ++idx)
    dst[idx] = data[idx];
}

template <typename T, std::size_t N, std::size_t... Idxes>
std::array<T, N> array_construct_aux(T const *src,
                                     std::index_sequence<Idxes...>) {
  return {src[Idxes]...};
}

template <typename T, std::size_t N>
std::array<T, N> array_construct(T const *src) {
  return array_construct_aux(src, std::make_index_sequence<N>{});
}

template <typename T, std::size_t N, typename Idx, std::size_t... I>
std::array<T, N> array_gather_aux(T const *src, Idx const *idxes,
                                  std::index_sequence<I...>) {
  return {src[idxes[I]]...};
}

template <typename T, std::size_t N, typename Idxes>
std::array<T, N> array_gather(T const *src, Idxes const &idxes) {
  lane_type_t<Idxes> idxes_[lane_size_v<Idxes>];
  vcl_store(idxes, idxes_);
  return array_gather_aux(src, idxes_, std::make_index_sequence<N>{});
}

template <typename T, std::size_t N, typename Idxes>
void array_scatter(std::array<T, N> const &data, T *dst, Idxes const &idxes) {
  lane_type_t<Idxes> idxes_[lane_size_v<Idxes>];
  vcl_store(idxes, idxes_);
  for (int idx = 0; idx < N; ++idx)
    dst[idxes_[idx]] = data[idx];
}

template <typename T, std::size_t N, typename Idx, std::size_t... I>
std::array<T, N> array_lookup_(T const *src, Idx const *idxes,
                               std::index_sequence<I...>) {
  return {src[idxes[I]]...};
}

template <typename T, std::size_t N, typename Idx>
std::array<T, N> array_lookup(T const *src, Idx const *idxes) {
  return array_lookup_(src, idxes, std::make_index_sequence<N>{});
}
} // namespace nitro::def