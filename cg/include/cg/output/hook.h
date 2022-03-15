#pragma once
#include <ioxx/ioxx.h>

namespace cg::out {
struct report_state {
  bool simul_first_time, traj_first_time;
  int *traj_idx, step_idx;
  std::filesystem::path out_dir;
  ioxx::xyaml::node for_simul, for_traj, for_step;

  void on_new_trajectory();
};

class hook {
public:
  virtual ~hook() = default;
  virtual void report_to(report_state &report) const = 0;
};
} // namespace cg::out