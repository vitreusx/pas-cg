#include <cg/pauli/update_pairs.h>
namespace cg::pauli {

void update_pairs::operator()() const {
  pairs->clear();

  auto cutoff = r_excl;

  for (int pair_idx = 0; pair_idx < nl->non_native.size(); ++pair_idx) {
    auto pair = nl->non_native[pair_idx];
    auto i1 = pair.i1(), i2 = pair.i2();

    auto r1 = r[i1], r2 = r[i2];
    auto dist = norm(simul_box->wrap(r1, r2));
    if (dist < cutoff + nl->orig_pad) {
      pairs->emplace_back(i1, i2);
    }
  }
}
} // namespace cg::pauli