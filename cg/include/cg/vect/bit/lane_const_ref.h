#pragma once
#include "../def/lane.h"
#include "lane_io.h"

namespace nitro::bit {
template <std::size_t N, std::size_t W>
class lane_const_ref {
public:
  explicit lane_const_ref(byte const *p) : p{p} {
    nitro::bit::load(_lane, p);
  }

  operator def::lane<bool, N, W>() const {
    return _lane;
  }

  template <typename U>
  auto &operator=(U const &val) const = delete;

private:
  byte const *p;
  def::lane<bool, N, W> _lane;
};
}