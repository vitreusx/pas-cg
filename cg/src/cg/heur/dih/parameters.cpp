#include <cg/heur/dih/parameters.h>
#include <cg/utils/ioxx_interop.h>
namespace cg::heur_dih {

void parameters::load(ioxx::xyaml::node const &p) {
  p["enabled"] >> enabled;
  auto coeffs_csv = p["coefficients"].as<ioxx::xyaml::table_file>();
  for (auto const &row : coeffs_csv->rows) {
    pair_coeffs pc;
    row["type2"] >> pc.type2;
    row["type3"] >> pc.type3;
    row["const"] >> pc.const_.assumed_("eps");
    row["sin"] >> pc.sin.assumed_("eps");
    row["cos"] >> pc.cos.assumed_("eps");
    row["sin2"] >> pc.sin2.assumed_("eps");
    row["cos2"] >> pc.cos2.assumed_("eps");
    row["sin_cos"] >> pc.sin_cos.assumed_("eps");

    coeffs[aa_heur_pair(pc.type2, pc.type3)] = pc;
  }
}
} // namespace cg::heur_dih