#include <cg/nl/parameters.h>

namespace cg::nl {

void parameters::load(ioxx::xyaml::node const &node) {
  std::string algorithm_str;
  node["algorithm"] >> algorithm_str;
  if (algorithm_str == "cell")
    algorithm = CELL;
  else
    algorithm = LEGACY;

  node["pad"] >> pad;

  if (node["cutoff"].Scalar() != "derived")
    node["cutoff"] >> cutoff;
}
} // namespace cg::nl