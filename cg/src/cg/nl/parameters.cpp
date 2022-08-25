#include <cg/nl/parameters.h>

namespace cg::nl {

void parameters::load(ioxx::xyaml::node const &node) {
  node["pad"] >> pad;

  if (node["cutoff"].Scalar() != "derived")
    node["cutoff"] >> cutoff;
}
} // namespace cg::nl