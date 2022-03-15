#pragma once
#include <cg/types/amp.h>
#include <ioxx/ioxx.h>

namespace cg::out {
struct report_data {
  bool simul_first_time, traj_first_time;
  int *traj_idx, step_idx;
  std::filesystem::path out_dir, traj_dir, step_dir;
  ioxx::xyaml::node for_simul, for_traj, for_step;
  bool report_stats, report_files;
  real last_stats_t, last_files_t;

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
    scalars.back()[key] = ioxx::convert<std::string>(value);
  }

  template <typename T>
  void add_step_scalar(std::string const &key, ioxx::xyaml::node &out,
                       T const &value) {
    out = value;
    add_step_scalar(key, value);
  }
};
} // namespace cg::out