#pragma once
#include <string>
#include <vector>

namespace cg::out {
class chain_id_seq_ {
public:
  chain_id_seq_();
  std::string const &operator[](int chain_idx) const;

private:
  std::vector<std::string> ids;
};
} // namespace cg::out