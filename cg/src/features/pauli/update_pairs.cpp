#include "features/pauli/update_pairs.h"
using namespace cg::pauli;

void update_pairs::operator()() const {
  pairs->clear();

  auto cutoff = r_excl;
  auto min_norm_inv = (real)1.0 / (cutoff + nl->orig_pad);

  for (int pair_idx = 0; pair_idx < nl->non_native.size(); ++pair_idx) {
    auto pair = nl->non_native[pair_idx];
    auto i1 = pair.i1(), i2 = pair.i2();

    auto r1 = r->at(i1), r2 = r->at(i2);
    if (norm_inv(box->r_uv(r1, r2)) > min_norm_inv) {
      pairs->emplace_back(i1, i2);
    }
  }
}