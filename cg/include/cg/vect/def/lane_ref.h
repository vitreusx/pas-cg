#pragma once
#include "lane.h"
#include "lane_io.h"

namespace nitro::def {
template <typename T, std::size_t N, std::size_t W>
class lane_ref {
public:
  explicit lane_ref(T *p) : p{p} {
    load(_lane, p);
  }

  operator lane<T, N, W>() const {
    return _lane;
  }

  template <typename U>
  auto &operator=(U const &val) const {
    _lane = val;
    store(_lane, p);
    return *this;
  }

private:
  T *p;
  mutable lane<T, N, W> _lane;
};
} // namespace nitro::def