#include "heur/dih/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::heur_dih;

void parameters::pair_coeffs::connect(ioxx::row_proxy &proxy) {
  proxy["type2"] & type2;
  proxy["type3"] & type3;
  proxy["const"] & const_.in_("eps");
  proxy["sin"] & sin.in_("eps");
  proxy["cos"] & cos.in_("eps");
  proxy["sin2"] & sin2.in_("eps");
  proxy["cos2"] & cos2.in_("eps");
  proxy["sin_cos"] & sin_cos.in_("eps");
}

void parameters::load(ioxx::xyaml::node const &p) {
  enabled = p["enabled"].as<bool>();
  auto coeffs_csv = p["coefficients"].as<ioxx::csv<pair_coeffs>>();
  for (auto const &row : coeffs_csv.rows) {
    coeffs[aa_heur_pair(row.type2, row.type3)] = row;
  }
}