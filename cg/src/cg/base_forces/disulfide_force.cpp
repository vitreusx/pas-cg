#include <cg/base_forces/disulfide.h>
#include <cg/utils/ioxx_interop.h>
namespace cg {

void disulfide_force::load(ioxx::xyaml::node const &node) {
  if (auto harm_node = node["harmonic"]; harm_node) {
    auto H1 = harm_node["H1"].as<quantity>();
    auto H2 = harm_node["H2"].as<quantity>();
    auto equil_dist = harm_node["equilibrium dist"].as<quantity>();
    force = harmonic(H1, H2, equil_dist);
  } else if (auto lj_node = node["lj"]; lj_node) {
    auto depth = lj_node["depth"].as<quantity>();
    auto r_min = lj_node["r_min"].as<quantity>();
    force = lj(depth, r_min);
  }
}
} // namespace cg