#include <cg/qa/legacy.h>

namespace cg::qa {
void legacy_impl::update_verlet_list() const {
  kqist[jq].clear();
  int kqactual = 0, kqont2 = kqist[jq2].size();

  for (int pair_idx = 0; pair_idx < nl->non_native.size(); ++pair_idx) {
    auto nl_pair = nl->non_native[pair_idx];
    auto i1 = nl_pair.i1(), i2 = nl_pair.i2();
    auto r1 = r[i1], r2 = r[i2];

    auto chain1 = chain_idx[i1], chain2 = chain_idx[i2];
    auto seq1 = seq_idx[i1], seq2 = seq_idx[i2];
    if (!include4 && chain1 == chain2 && abs(seq1 - seq2) == 4)
      continue;

    auto dist = norm(simul_box->wrap(r1, r2));
    if (dist < cutoff + nl->orig_pad) {
      bool found = false;
      while (kqactual < kqont2) {
        auto kq1 = kqist[jq2][kqactual][0], kq2 = kqist[jq2][kqactual][1];
        if (kq1 < i1 || (kq1 == i1 && kq2 < i2)) {
          ++kqactual;
        } else if (kq1 == i1 && kq2 == i2) {
          auto const &prev = kqist[jq2][kqactual];
          qa_data_t cur = {i1, i2, prev[2], prev[3]};
          kqist[jq].push_back(cur);
          found = true;
          ++kqactual;
          break;
        } else {
          break;
        }
      }

      if (!found) {
        auto lsamechain = chain_idx[i1] == chain_idx[i2];
        qa_data_t cur = {i1, i2, 0, lsamechain ? 1 : -1};
        kqist[jq].push_back(cur);

        if (charge[(uint8_t)atype[i1]] != 0 && charge[(uint8_t)atype[i2]] != 0) {

        }
      }
    }
  }
}

void legacy_impl::evalcpot() const {}
} // namespace cg::qa