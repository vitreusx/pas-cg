#include "qa/update_free_pairs.h"
using namespace cg::qa;

void update_free_pairs::operator()() const {
  pairs->clear();

  auto cutoff = max_formation_min_dist;
  auto min_norm_inv = (real)1.0 / (cutoff + nl->orig_pad);

  for (int pair_idx = 0; pair_idx < nl->non_native.size(); ++pair_idx) {
    auto nl_pair = nl->non_native[pair_idx];
    auto i1 = nl_pair.i1(), i2 = nl_pair.i2();
    auto r1 = r[i1], r2 = r[i2];

    if (norm_inv(simul_box->r_uv(r1, r2)) > min_norm_inv) {
      pairs->emplace(i1, i2);
    }
  }
}