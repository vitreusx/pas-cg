#include <cg/amino/aa_data.h>

#include <cg/utils/quantity.h>
#include <numeric>

namespace cg {
void contact_limits::link(ioxx::xyaml::proxy &proxy) {
  proxy["back"] & back;
  proxy["side (all)"] & side_all;
  proxy["side (hydrophobic)"] & side_hydrophobic;
  proxy["side (polar)"] & side_polar;
}

static bool close(double x, double y) {
  return abs(x - y) / abs(x) < 1e-5;
}

static double correct(std::string const &name, double radius) {
  if (name[0] == 'N')
    return quantity(1.64, "A");
  else if (name[0] == 'S')
    return quantity(1.77, "A");
  else if (name[0] == 'O' && !close(radius, quantity(1.42, "A")) &&
           !close(radius, quantity(1.46, "A")))
    return quantity(1.46, "A");
  else if (name[0] == 'C' && !close(radius, quantity(1.88, "A")) &&
           !close(radius, quantity(1.76, "A")) &&
           !close(radius, quantity(1.61, "A")))
    return quantity(1.88, "A");
  else
    return radius;
}

void amino_acid_data::load(ioxx::xyaml::node const &node) {
  using namespace ioxx::xyaml;

  subnode sub;
  node >> sub;

  for (auto const &entry : sub["default atom data"]) {
    auto name = entry.first.as<std::string>();
    auto &atom_data = def_atoms[name];
    atom_data.name = name;
    atom_data.radius = quantity(entry.second.as<std::string>()).assumed_("A");
    atom_data.radius = correct(atom_data.name, atom_data.radius);
  }

  for (auto const &entry : sub["amino acids"]) {
    auto name = entry.first.as<std::string>();
    auto data_node = ioxx::xyaml::node(entry.second);

    aa_data &cur_data = data[amino_acid(name)];
    cur_data.mass = data_node["mass"].as<cg::quantity>().assumed("amu");
    cur_data.radius = data_node["radius"].as<cg::quantity>().assumed("A");

    if (auto alt = data_node["alt atom data"]; alt) {
      for (auto const &alt_entry : alt) {
        auto atom_name = alt_entry.first.Scalar();
        auto &override = cur_data.overrides[atom_name];
        override.name = atom_name;
        override.radius = quantity(alt_entry.second.Scalar()).assumed_("A");
        override.radius = correct(override.name, override.radius);
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

atom_data const &amino_acid_data::for_atom(const amino_acid &res,
                                           const std::string &name) const {
  for (int n = (int)name.size(); n > 0; --n) {
    auto prefix = name.substr(0, n);
    for (auto *map : {&data.at(res).overrides, &def_atoms}) {
      if (auto iter = map->find(prefix); iter != map->end())
        return iter->second;
    }
  }

  std::stringstream error_ss;
  error_ss << "no data found for an atom \"" << name << "\" for residue \""
           << res.name() << "\"";
  throw std::runtime_error(error_ss.str());
}
} // namespace cg