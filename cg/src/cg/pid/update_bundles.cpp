#include <cg/pid/update_bundles.h>
namespace cg::pid {

void update_bundles::operator()() const {
  bundles->clear();
  end_bundles.clear();

  auto n = r.size();
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
      auto i1p = i1 - 1, i1n = i1 + 1, i2p = i2 - 1, i2n = i2 + 1;
      auto val_i1p = (i1p >= 0), val_i1n = (i1n < n), val_i2p = (i2p >= 0),
           val_i2n = (i2n < n);
      auto terminal = !(val_i1p && val_i1n && val_i2p && val_i2n);
      auto *bundles_ = terminal ? &end_bundles : bundles;

      auto atype1 = atype[i1], atype2 = atype[i2];

      int16_t type = (int16_t)(uint8_t)atype1 * (int16_t)amino_acid::NUM_TYPES +
                     (int16_t)(uint8_t)atype2;
      bundles_->emplace_back(i1, i2, nl_pair.orig_dist(), type);
      //      nl_pair.taken() = true;
    }
  }

  std::sort(
      bundles->begin(), bundles->end(),
      [](auto const &x, auto const &y) -> auto{
        return x.orig_dist() < y.orig_dist();
      });

  *fast_iter_end = bundles->size();
  for (auto bundle : end_bundles)
    bundles->push_back(bundle);
}
} // namespace cg::pid