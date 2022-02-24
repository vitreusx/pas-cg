#include "amino/compiled.h"
using namespace cg;

compiled_aa_data::compiled_aa_data(amino_acid_data const &data) {
  mass = charge = nitro::vector<real>(amino_acid::NUM_TYPES);
  ptype = nitro::vector<polarization_type>(amino_acid::NUM_TYPES);
  base_sync = nitro::vector<sync_data>(amino_acid::NUM_TYPES);

  for (auto const &[code, data_for_amino] : data.data) {
    auto aa_idx = (uint8_t)code;
    mass[aa_idx] = data_for_amino.mass;
    charge[aa_idx] = data_for_amino.charge;
    ptype[aa_idx] = data_for_amino.polarization;

    auto const &lim = data_for_amino.limits;
    base_sync[aa_idx] =
        sync_data(lim.back, lim.side_all, lim.side_polar, lim.side_hydrophobic);
  }
}