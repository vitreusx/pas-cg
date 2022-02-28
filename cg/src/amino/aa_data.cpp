#include "amino/aa_data.h"
#include "utils/ioxx_interop.h"
#include "utils/quantity.h"
#include <numeric>
using namespace cg;

void contact_limits::link(ioxx::xyaml::proxy &proxy) {
  proxy["back"] & back;
  proxy["side (all)"] & side_all;
  proxy["side (hydrophobic)"] & side_hydrophobic;
  proxy["side (polar)"] & side_polar;
}

void amino_acid_data::load(ioxx::xyaml::node const &node) {
  using namespace ioxx::xyaml;

  subnode sub;
  node >> sub;

  std::unordered_map<std::string, double> def_atom_radii;
  for (auto const &entry : sub["default atom radii"]) {
    auto name = entry.first.as<std::string>();
    auto radius = quantity(entry.second.as<std::string>()).assumed_("A");
    def_atom_radii[name] = radius;
  }

  for (auto const &entry : sub["amino acids"]) {
    auto name = entry.first.as<std::string>();
    auto data_node = ioxx::xyaml::node(entry.second);

    aa_data &cur_data = data[amino_acid(name)];
    data_node["mass"] >> cur_data.mass.assumed_("amu");
    data_node["radius"] >> cur_data.radius.assumed_("A");

    auto atom_radii = def_atom_radii;
    if (auto alt = data_node["alt atom radii"]; alt) {
      for (auto const &alt_entry : alt) {
        auto atom_name = alt_entry.first.Scalar();
        auto alt_radius = quantity(alt_entry.second.Scalar()).assumed_("A");
        atom_radii[atom_name] = alt_radius;
      }
    }

    for (auto const &back_atom : {"N", "CA", "C", "O", "OXT"}) {
      atom_data &atom_ = cur_data.atoms[back_atom];
      atom_.name = back_atom;
      atom_.radius = atom_radii[back_atom];
      atom_.backbone = true;
    }

    if (auto side_atoms_node = data_node["side"]; side_atoms_node) {
      for (auto const &side_atom_node : side_atoms_node) {
        auto side_atom = side_atom_node.as<std::string>();
        atom_data &atom_ = cur_data.atoms[side_atom];
        atom_.name = side_atom;
        atom_.radius = atom_radii[side_atom];
        atom_.backbone = false;
      }
    }

    cur_data.polarization = polarization_type::NONE;
    if (auto polarization_node = data_node["polarization"]; polarization_node) {
      auto ptype = polarization_node.as<std::string>();
      if (ptype == "polar") {
        cur_data.polarization = polarization_type::POLAR;
      } else if (ptype == "hydrophobic") {
        cur_data.polarization = polarization_type::HYDROPHOBIC;
      }
    }

    cur_data.charge = 0.0;
    if (auto charge_node = data_node["charge"]; charge_node) {
      auto charge = quantity(charge_node.as<std::string>()).assumed_("e");
      cur_data.charge = charge;
    }

    data_node["contact limits"] >> cur_data.limits;
  }

  auto avg_mass = std::accumulate(
      data.begin(), data.end(), 0.0,
      [](auto const &sum, auto const &entry) -> auto {
        auto const &[name, res_data] = entry;
        return sum + res_data.mass;
      });
  avg_mass /= (double)data.size();

  for (auto &[name, res_data] : data) {
    res_data.mass = quantity(res_data.mass / avg_mass, "f77mass");
  }
}

aa_data const &amino_acid_data::operator[](const amino_acid &aa) const {
  return data.at(aa);
}