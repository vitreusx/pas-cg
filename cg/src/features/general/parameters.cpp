#include "features/general/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::gen;

void parameters::connect(ioxx::xyaml_proxy &p) {
  p["total time"] & total_time;
  p["seed"] & seed;
  p["num of threads"] & num_of_threads;
}