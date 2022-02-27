#include "nl/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::nl;

void parameters::load(ioxx::xyaml::node const &node) {
  std::string algorithm_str;
  node["algorithm"] >> algorithm_str;
  if (algorithm_str == "cell")
    algorithm = CELL;
  else
    algorithm = LEGACY;
  node["pad factor"] >> pad_factor;
}