#include "heur/ang/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::heur_ang;

void parameters::pair_coeffs::connect(ioxx::row_proxy &proxy) {
  proxy["type1"] & type1;
  proxy["type2"] & type2;
  for (int deg = 0; deg <= 6; ++deg)
    proxy["a" + std::to_string(deg)] & poly[deg].assumed_("eps");
}

void parameters::load(ioxx::xyaml::node const &p) {
  enabled = p["enabled"].as<bool>();
  auto coeffs_csv = p["coefficients"].as<ioxx::csv<pair_coeffs>>();
  for (auto const &row : coeffs_csv.rows) {
    coeffs[aa_heur_pair(row.type1, row.type2)] = row;
  }
}