#include <cg/general/parameters.h>
#include <cg/utils/ioxx_interop.h>
namespace cg::gen {

void parameters::load(ioxx::xyaml::node const &p) {
  p["total time"] >> total_time;
  p["equil time"] >> equil_time;
  p["seed"] >> seed;
  p["num of threads"] >> num_of_threads;
  p["num of trajectories"] >> num_of_traj;
  p["disable all forces"] >> disable_all;

  auto dp = p["debug mode"];
  dp["enabled"] >> debug_mode.enabled;
  dp["floating point exceptions"] >> debug_mode.fp_exceptions;
  dp["determinism"] >> debug_mode.determinism;
  debug_mode.fp_exceptions &= debug_mode.enabled;
  debug_mode.determinism &= debug_mode.enabled;
}
} // namespace cg::gen