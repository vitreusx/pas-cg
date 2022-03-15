#pragma once
#include <ioxx/ioxx.h>

namespace cg::out {
struct report_state {
  bool simul_first_time, traj_first_time;
  int *traj_idx, step_idx;
  std::filesystem::path out_dir, traj_dir, step_dir;
  ioxx::xyaml::node for_simul, for_traj, for_step;

  std::vector<std::string> columns;
  std::vector<std::unordered_map<std::string, std::string>> scalars;

  void simul_init();
  void traj_init();
  void step_init();
  void step_finish();

  template <typename T>
  void add_step_scalar(std::string const &key, T const &value) {
    if (std::find(columns.begin(), columns.end(), key) == columns.end())
      columns.push_back(key);
    scalars.back()[key] = value;
  }
};
} // namespace cg::out