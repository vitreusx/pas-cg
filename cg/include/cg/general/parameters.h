#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::gen {
struct parameters {
  quantity total_time, equil_time;
  int seed, num_of_threads, num_of_traj;
  bool debug_mode;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::gen