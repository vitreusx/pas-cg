#include <cg/dh/parameters.h>

namespace cg::dh {
void const_dh_parameters::load(const ioxx::xyaml::node &node) {
  node["permittivity"] >> permittivity;
}

void rel_dh_parameters::load(const ioxx::xyaml::node &node) {
  node["permittivity factor"] >> perm_factor;
}

void parameters::load(const ioxx::xyaml::node &node) {
  node["enabled"] >> enabled;
  node["variant"] >> variant;
  node["screening distance"] >> screening_dist;
  node["cutoff distance"] >> cutoff;
  node["constant variant params"] >> const_dh;
  node["relative variant params"] >> rel_dh;
}
} // namespace cg::dh