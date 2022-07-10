#include <cg/base_forces/spec.h>

namespace cg {
void force_spec::load(const ioxx::xyaml::node &node) {
  if (node["variant"])
    node["variant"] >> variant;

  node["harmonic params"] >> harmonic;
  node["lj params"] >> lj;
  node["shifted lj params"] >> shifted_lj;
  node["sink lj params"] >> sink_lj;
}

void ss_force_spec::load(const ioxx::xyaml::node &node) {
  if (node["variant"])
    node["variant"] >> variant;

  node["lj params"] >> lj;
  node["shifted lj params"] >> shifted_lj;
  node["sink lj params"] >> sink_lj;
}
}