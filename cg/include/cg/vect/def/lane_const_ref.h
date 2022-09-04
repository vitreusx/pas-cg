#pragma once
#include "lane.h"
#include "lane_io.h"

namespace nitro::def {
template <typename T, std::size_t N, std::size_t W>
class lane_const_ref {
public:
  explicit lane_const_ref(T const *p) {
    load(_lane, p);
  }

  template <typename U>
  auto &operator=(U const &val) const = delete;

  operator lane<T, N, W>() const {
    return _lane;
  }

private:
  lane<T, N, W> _lane;
};
} // namespace nitro::def