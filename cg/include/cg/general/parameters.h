#pragma once
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg::gen {
struct parameters {
  enum class prog_mode {
    perform_simulation,
    check_determinism
  };

  struct debug_mode_t {
    bool enabled, fp_exceptions, disable_all;
  };

  prog_mode mode;
  quantity total_time, equil_time, repulsive_cutoff;
  std::optional<quantity> fixed_cutoff;
  int seed, num_of_threads, num_of_traj;
  bool pbc_x, pbc_y, pbc_z;

  debug_mode_t debug_mode;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::gen