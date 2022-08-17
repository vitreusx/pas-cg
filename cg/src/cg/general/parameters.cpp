#include <cg/general/parameters.h>

namespace cg::gen {
void parameters::load(ioxx::xyaml::node const &p) {
  p["total time"] >> total_time;
  p["equil time"] >> equil_time;
  p["seed"] >> seed;
  p["num of threads"] >> num_of_threads;
  p["num of trajectories"] >> num_of_traj;
  p["repulsive cutoff"] >> repulsive_cutoff;
  p["counting factor"] >> counting_factor;

  auto dp = p["debug mode"];
  dp["floating point exceptions"] >> fp_exceptions;
  dp["dump data for every step"] >> dump_data;

  auto mode_str = p["mode"].as<std::string>();
  if (mode_str == "perform simulation")
    mode = prog_mode::perform_simulation;
  else if (mode_str == "check determinism")
    mode = prog_mode::check_determinism;
}
} // namespace cg::gen