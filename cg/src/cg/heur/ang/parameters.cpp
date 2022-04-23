#include <cg/heur/ang/parameters.h>
#include <cg/utils/ioxx_interop.h>
namespace cg::heur_ang {

void parameters::load(ioxx::xyaml::node const &p) {
  enabled = p["enabled"].as<bool>();
  auto coeffs_csv = p["coefficients"].as<ioxx::xyaml::table_file>();
  for (auto const &row : coeffs_csv->rows) {
    pair_coeffs pc;
    row["type1"] >> pc.type1;
    row["type2"] >> pc.type2;
    for (int deg = 0; deg <= 6; ++deg) {
      auto unit = "eps/rad**" + std::to_string(deg);
      row["a" + std::to_string(deg)] >> pc.poly[deg].assumed_(unit);
    }

    coeffs[aa_heur_pair(pc.type1, pc.type2)] = pc;
  }
}
} // namespace cg::heur_ang