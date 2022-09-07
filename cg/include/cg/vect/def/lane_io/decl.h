#pragma once

namespace nitro::def {
template <typename Lane, typename T>
void load(Lane &data, T const *src);

template <typename Lane, typename T>
void store(Lane const &data, T *dst);

template <typename T>
struct selector {};

template <typename Lane, typename T>
Lane construct(T const *src) {
  return _construct(selector<Lane>{}, src);
}

template <typename Lane, typename T, typename Idxes>
Lane gather(T const *src, Idxes const &idxes) {
  return _gather(selector<Lane>{}, src, idxes);
}

template <typename Lane, typename T, typename Idxes>
void scatter(Lane const &data, T *dst, Idxes const &idxes);
} // namespace nitro::def