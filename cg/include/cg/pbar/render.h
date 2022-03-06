#pragma once
#include <cg/types/amp.h>
#include <chrono>

namespace cg::pbar {
class render {
public:
  int width;
  real total_time, period_s;

public:
  real const *V, *t;
  bool *is_first;
  using time_point_t = decltype(std::chrono::high_resolution_clock::now());
  time_point_t *start_clock, *last_clock;

public:
  void operator()() const;
};
} // namespace cg::pbar