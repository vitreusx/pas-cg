#include <cg/angles/parameters.h>

namespace cg::angles {
void heur_ang_parameters::load(const ioxx::xyaml::node &p) {
  enabled = p["enabled"].as<std::optional<bool>>().value_or(enabled);

  auto coeffs_csv = p["coefficients"].as<ioxx::xyaml::table_file>();
  for (auto const &row : coeffs_csv->rows) {
    pair_coeffs pc;
    row["type1"] >> pc.type1;
    row["type2"] >> pc.type2;
    for (int deg = 0; deg <= 6; ++deg) {
      auto unit = "eps/rad**" + std::to_string(deg);
      row["a" + std::to_string(deg)] >> pc.poly[deg].assumed(unit);
    }

    coeffs[aa_heur_pair(pc.type1, pc.type2)] = pc;
  }
}

void heur_dih_parameters::load(const ioxx::xyaml::node &p) {
  enabled = p["enabled"].as<std::optional<bool>>().value_or(enabled);

  auto coeffs_csv = p["coefficients"].as<ioxx::xyaml::table_file>();
  for (auto const &row : coeffs_csv->rows) {
    pair_coeffs pc;
    row["type2"] >> pc.type2;
    row["type3"] >> pc.type3;
    row["const"] >> pc.const_.assumed("eps");
    row["sin"] >> pc.sin.assumed("eps");
    row["cos"] >> pc.cos.assumed("eps");
    row["sin2"] >> pc.sin2.assumed("eps");
    row["cos2"] >> pc.cos2.assumed("eps");
    row["sin_cos"] >> pc.sin_cos.assumed("eps");

    coeffs[aa_heur_pair(pc.type2, pc.type3)] = pc;
  }
}

void nat_ang_parameters::load(const ioxx::xyaml::node &p) {
  enabled = p["enabled"].as<std::optional<bool>>().value_or(enabled);

  CBA = p["CBA"].as<quantity>();
}

void nat_dih_parameters::complex_variant_parameters::load(
    const ioxx::xyaml::node &p) {
  CDA = p["CDA"].as<quantity>();
  CDB = p["CDB"].as<quantity>();
}

void nat_dih_parameters::simple_variant_parameters::load(
    const ioxx::xyaml::node &p) {
  CDH = p["CDH"].as<quantity>();
}

void nat_dih_parameters::load(const ioxx::xyaml::node &node) {
  enabled = node["enabled"].as<std::optional<bool>>().value_or(enabled);
  node["variant"] >> variant;
  node["complex variant params"] >> complex;
  node["simple variant params"] >> simple;
}

void parameters::load(const ioxx::xyaml::node &node) {
  node["all enabled"] >> all_enabled;
  node["dihedral potentials enabled"] >> dih_enabled;

  heur_ang.enabled = all_enabled;
  nat_ang.enabled = all_enabled;
  heur_dih.enabled = all_enabled && dih_enabled;
  nat_dih.enabled = all_enabled && dih_enabled;

  node["heurestic bond angles params"] >> heur_ang;
  node["native bond angles params"] >> nat_ang;
  node["heurestic dihedral angles params"] >> heur_dih;
  node["native dihedral angles params"] >> nat_dih;
}
} // namespace cg::angles