#pragma once
#include <cg/types/amp.h>

namespace cg {
template <typename T> void sanitize(T &value) {
  if (!cg::isfinite(value))
    value = (T)0;
  else
    value = cg::clamp(value, (T)-1e3, (T)1e3);
}

template <typename T> void sanitize(cg::vec3<T> &v) {
  sanitize(v.x());
  sanitize(v.y());
  sanitize(v.z());
}
} // namespace cg