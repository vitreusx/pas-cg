#include <cg/output/chain_id_seq.h>
#include <cg/utils/text.h>

namespace cg::out {
chain_id_seq_::chain_id_seq_() {
  ids.clear();

  for (char c = 'A'; c <= 'Z'; ++c)
    ids.push_back(format(" %c", c));

  for (char c = 'a'; c <= 'z'; ++c)
    ids.push_back(format(" %c", c));

  for (int idx = 0; idx < 100; ++idx)
    ids.push_back(format("%2d", idx));

  ids.emplace_back("**");
}

std::string const &chain_id_seq_::operator[](int chain_idx) const {
  return ids[chain_idx % ids.size()];
}
} // namespace cg::out