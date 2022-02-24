#pragma once
#include <cg/amino/amino_acid.h>

namespace cg {
class aa_heur_pair {
public:
  aa_heur_pair() = default;
  explicit aa_heur_pair(char type1, char type2);
  explicit aa_heur_pair(amino_acid const &a1, amino_acid const &a2);
  operator uint8_t() const;

private:
  uint8_t value;
};
} // namespace cg

namespace std {
template <> struct hash<cg::aa_heur_pair> {
  size_t operator()(cg::aa_heur_pair const &pair) const;
};
} // namespace std