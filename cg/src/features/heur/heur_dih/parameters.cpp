#include "features/heur/heur_dih/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::heur_dih;

void parameters::pair_coeffs::connect(ioxx::row_proxy &proxy) {
  proxy["type2"] & type2;
  proxy["type3"] & type3;
  proxy["sin"] & sin;
  proxy["cos"] & cos;
  proxy["sin2"] & sin2;
  proxy["cos2"] & cos2;
  proxy["sin_cos"] & sin_cos;
}

void parameters::connect(ioxx::xyaml_proxy &p) {
  enabled = p["enabled"].as<bool>();
  auto coeffs_csv = p["coefficients"].as<ioxx::csv<pair_coeffs>>();
  for (auto const &row : coeffs_csv.rows) {
    coeffs[aa_heur_pair(row.type2, row.type3)] = row;
  }
}