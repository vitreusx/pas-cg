#pragma once
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>

namespace cg::gen {
struct parameters {
  struct debug_mode_t {
    bool enabled, fp_exceptions, determinism;
  };

  quantity total_time, equil_time;
  int seed, num_of_threads, num_of_traj;
  bool disable_all;
  
  debug_mode_t debug_mode;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::gen