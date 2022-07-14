#include <cg/base_forces/harmonic.h>


namespace cg {
void harmonic_specs::load(const ioxx::xyaml::node &node) {
  node["H1"] >> H1;
  node["H2"] >> H2;
  node["nat_r"] >> nat_r;
}

harmonic_specs::operator harmonic() const {
  return harmonic(H1.value(), H2.value(), nat_r.value());
}
} // namespace cg