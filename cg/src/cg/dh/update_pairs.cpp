#include <cg/dh/update_pairs.h>
namespace cg::dh {

void update_pairs::operator()() const {
  pairs->clear();

  auto min_norm_inv = (real)1.0 / (cutoff + nl->orig_pad);

  for (int pair_idx = 0; pair_idx < nl->non_native.size(); ++pair_idx) {
    auto pair = nl->non_native.at(pair_idx);
    auto i1 = pair.i1(), i2 = pair.i2();
    auto r1 = r[i1], r2 = r[i2];

    if (norm_inv(simul_box->wrap(r1, r2)) > min_norm_inv) {
      auto q1_x_q2 = q[(uint8_t)atype[i1]] * q[(uint8_t)atype[i2]];
      if (q1_x_q2 != 0.0) {
        pairs->emplace_back(i1, i2, q1_x_q2);
      }
    }
  }
}
} // namespace cg::dh