#include <cg/input/morph_into_saw.h>
#include <cg/utils/ioxx_interop.h>
namespace cg::input {

void morph_into_saw_t::link(ioxx::xyaml::proxy &p) {
  p["perform"] & perform;
  p["residue density"] & residue_density;
  p["intersection at"] & intersection_at;
  p["num of retries"] & num_of_retries;
  p["sample initial positions from the box"] & sample_from_box;
  p["PBC also during SAW"] & with_pbc;

  auto bond_dist_str = p["bond distance"].as<std::optional<std::string>>();
  if (bond_dist_str.has_value() && bond_dist_str.value() != "average") {
    bond_distance = quantity(bond_dist_str.value());
  } else {
    bond_distance = std::nullopt;
  }
}
} // namespace cg::input