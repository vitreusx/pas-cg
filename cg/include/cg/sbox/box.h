#pragma once
#include <cg/types/vec3.h>
#include <cg/utils/math.h>

namespace cg::sbox {
template <typename U> class box {
public:
  template <typename E> inline void extend(vec3_expr<E> const &v) {
    if (min.x() > max.x()) {
      min = max = v;
    } else {
      min.x() = cg::min(min.x(), v.x());
      min.y() = cg::min(min.y(), v.y());
      min.z() = cg::min(min.z(), v.z());
      max.x() = cg::max(max.x(), v.x());
      max.y() = cg::max(max.y(), v.y());
      max.z() = cg::max(max.z(), v.z());
    }
  }

  auto extent() const {
    return max - min;
  }

  auto center() const {
    return (max + min) / (U)2.0;
  }

public:
  vec3<U> min = vec3<U>(1.0, 1.0, 1.0);
  vec3<U> max = vec3<U>(-1.0, -1.0, -1.0);
};
} // namespace cg::sbox