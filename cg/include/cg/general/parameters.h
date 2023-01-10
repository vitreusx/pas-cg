#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::gen {
struct parameters {
  enum class prog_mode {
    perform_simulation,
    check_determinism
  };

  prog_mode mode;
  quantity total_time, equil_time, repulsive_cutoff, counting_factor;
  int seed, num_of_threads, num_of_traj, dump_data_period;
  bool fp_exceptions, dump_data, use_gpu;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::gen