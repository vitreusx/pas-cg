#include <cg/general/parameters.h>
#include <cg/utils/ioxx_interop.h>
namespace cg::gen {

void parameters::load(ioxx::xyaml::node const &p) {
  p["total time"] >> total_time;
  p["equil time"] >> equil_time;
  p["seed"] >> seed;
  p["num of threads"] >> num_of_threads;
  p["num of trajectories"] >> num_of_traj;

  auto mode_str = p["mode"].as<std::string>();
  if (mode_str == "perform simulation")
    mode = prog_mode::perform_simulation;
  else if (mode_str == "check determinism")
    mode = prog_mode::check_determinism;
  else
    throw std::runtime_error("Unknown value for [general.mode]");

  auto dp = p["debug mode"];
  dp["enabled"] >> debug_mode.enabled;
  dp["floating point exceptions"] >> debug_mode.fp_exceptions;
  dp["disable all forces"] >> debug_mode.disable_all;
  debug_mode.fp_exceptions &= debug_mode.enabled;
  debug_mode.disable_all &= debug_mode.enabled;
}
} // namespace cg::gen