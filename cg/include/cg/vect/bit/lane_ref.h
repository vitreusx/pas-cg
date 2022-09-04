#pragma once
#include "../def/lane.h"
#include "lane_io.h"

namespace nitro::bit {
template <std::size_t N, std::size_t W>
class lane_ref {
public:
  explicit lane_ref(byte *p) : p{p} {
    load(_lane, p);
  }

  operator def::lane<bool, N, W>() const {
    return _lane;
  }

  template <typename U>
  auto &operator=(U const &val) const {
    _lane = val;
    store(_lane, p);
    return *this;
  }

private:
  byte *p;
  mutable def::lane<bool, N, W> _lane;
};
} // namespace nitro::bit