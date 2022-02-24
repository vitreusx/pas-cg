#pragma once
#include <cg/types/amp.h>
#include <chrono>

namespace cg::pbar {
class render {
public:
  int width;
  real total_time, period;

public:
  real *last_t;
  real const *V, *t;
  bool *is_first;
  using time_point_t = decltype(std::chrono::high_resolution_clock::now());
  time_point_t *start_wall_time;

public:
  void operator()() const;
};
} // namespace cg::pbar