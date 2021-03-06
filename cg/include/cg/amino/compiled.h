#pragma once
#include "aa_data.h"
#include "sync_data.h"
#include <cg/types/amp.h>

namespace cg {
class compiled_aa_data {
public:
  vect::vector<real> mass, charge;
  vect::vector<polarization_type> ptype;
  vect::vector<sync_data> base_sync;

public:
  compiled_aa_data() = default;
  explicit compiled_aa_data(amino_acid_data const &data);
};
} // namespace cg