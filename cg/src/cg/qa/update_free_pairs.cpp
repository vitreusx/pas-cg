#include <cg/qa/update_free_pairs.h>
namespace cg::qa {

void update_free_pairs::operator()() const {
  pairs->clear();

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
      pairs->emplace(i1, i2, nl_pair.orig_dist());
      nl_pair.taken() = true;
    }
  }
}
} // namespace cg::qa