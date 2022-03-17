#pragma once
#include <cg/types/amp.h>

namespace cg {
template <typename T> void sanitize(T &value, T const &max) {
  if (!cg::isfinite(value))
    value = (T)0;
  else
    value = cg::clamp(value, -max, max);
}

template <typename T> void sanitize(cg::vec3<T> &v, T const &max) {
  sanitize(v.x(), max);
  sanitize(v.y(), max);
  sanitize(v.z(), max);
}
} // namespace cg