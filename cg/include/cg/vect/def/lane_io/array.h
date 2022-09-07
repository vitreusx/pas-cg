#pragma once
#include "../lane/array.h"
#include "decl.h"

namespace nitro::def {

template <typename T, std::size_t N>
void load(std::array<T, N> &data, T const *src) {
  for (int idx = 0; idx < N; ++idx)
    data[idx] = src[idx];
}

template <typename T, std::size_t N>
void load(std::array<T, N> const &data, T *dst) {
  for (int idx = 0; idx < N; ++idx)
    dst[idx] = data[idx];
}

template <typename T, std::size_t N, std::size_t... Idxes>
std::array<T, N> _construct_aux(selector<std::array<T, N>>, T const *src,
                                std::index_sequence<Idxes...>) {
  return {src[Idxes]...};
}

template <typename T, std::size_t N>
std::array<T, N> _construct(selector<std::array<T, N>>, T const *src) {
  return _construct_aux(selector<std::array<T, N>>{}, src,
                        std::make_index_sequence<N>{});
}

template <typename T, std::size_t N, typename Idx, std::size_t... I>
std::array<T, N> _gather_aux(selector<std::array<T, N>>, T const *src,
                             Idx const *idxes, std::index_sequence<I...>) {
  return {src[idxes[I]]...};
}

template <typename T, std::size_t N, typename Idxes>
std::array<T, N> _gather(selector<std::array<T, N>>, T const *src,
                         Idxes const &idxes) {
  lane_type_t<Idxes> idxes_[lane_size_v<Idxes>];
  store(idxes, idxes_);
  return _gather_aux(selector<std::array<T, N>>{}, src, idxes_,
                     std::make_index_sequence<N>{});
}

template <typename T, std::size_t N, typename Idxes>
void scatter(std::array<T, N> const &data, T *dst, Idxes const &idxes) {
  lane_type_t<Idxes> idxes_[lane_size_v<Idxes>];
  store(idxes, idxes_);
  for (int idx = 0; idx < N; ++idx)
    dst[idxes_[idx]] = data[idx];
}
} // namespace nitro::def