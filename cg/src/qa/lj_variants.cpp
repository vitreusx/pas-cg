#include "qa/lj_variants.h"
namespace cg::qa {

lj_variants::lj_variants(const cg::lj_variants &variants) {
  this->variants.resize(contact_type::NUM_TYPES);

  this->variants[(int16_t)contact_type::BACK_BACK()] = sink_lj(variants.bb);
  this->variants[(int16_t)contact_type::BACK_SIDE()] = sink_lj(variants.bs);
  this->variants[(int16_t)contact_type::SIDE_BACK()] = sink_lj(variants.sb);

  for (auto const &aa1 : amino_acid::all()) {
    for (auto const &aa2 : amino_acid::all()) {
      this->variants[(int16_t)contact_type::SIDE_SIDE(aa1, aa2)] =
          variants.ss.at({aa1, aa2});
    }
  }
}
} // namespace cg::qa