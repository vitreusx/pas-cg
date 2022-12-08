#include <cg/pid/update_bundles.h>
namespace cg::pid {

void update_bundles::operator()() const {
  bundles->clear();
  end_bundles.clear();

  for (decltype(auto) nl_pair : nl->pairs) {
    if (nl_pair.taken())
      continue;

    auto i1 = nl_pair.i1(), i2 = nl_pair.i2();
    auto r1 = r[i1], r2 = r[i2];

    auto chain1 = chain_idx[i1], chain2 = chain_idx[i2];
    auto seq1 = seq_idx[i1], seq2 = seq_idx[i2];
    if (!include4 && chain1 == chain2 && abs(seq1 - seq2) == 4)
      continue;

    auto dist = norm(simul_box->wrap<vec3r>(r1, r2));
    if (dist < cutoff + nl->orig_pad) {
      auto prev1 = prev[i1], next1 = next[i1];
      auto prev2 = prev[i2], next2 = next[i2];

      auto terminal = prev1 < 0 || next1 < 0 || prev2 < 0 || next2 < 0;
      auto *bundles_ = terminal ? &end_bundles : bundles;
      //      auto *bundles_ = bundles;

      auto atype1 = atype[i1], atype2 = atype[i2];

      int16_t type = (int16_t)(uint8_t)atype1 * (int16_t)amino_acid::NUM_TYPES +
                     (int16_t)(uint8_t)atype2;
      bundles_->emplace_back(i1, i2, nl_pair.orig_dist(), type);
      //      nl_pair.taken() = true;
    }
  }

  *fast_iter_end = bundles->size();
  for (auto bundle : end_bundles)
    bundles->push_back(bundle);
}
} // namespace cg::pid