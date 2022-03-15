#include "general/parameters.h"
#include "utils/ioxx_interop.h"
namespace cg::gen {

void parameters::load(ioxx::xyaml::node const &p) {
  p["total time"] >> total_time;
  p["equil time"] >> equil_time;
  p["seed"] >> seed;
  p["num of threads"] >> num_of_threads;
  p["debug mode"] >> debug_mode;
  p["num of trajectories"] >> num_of_traj;
}
} // namespace cg::gen