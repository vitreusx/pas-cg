#include "features/heur/aa_heur_pair.h"
using namespace cg;

size_t std::hash<aa_heur_pair>::operator()(const cg::aa_heur_pair &pair) const {
  return std::hash<uint8_t>()(static_cast<uint8_t>(pair));
}

aa_heur_pair::aa_heur_pair(char type1, char type2) {
  auto ord1 = (type1 == 'G' ? 0u : (type1 == 'P' ? 1u : 2u));
  auto ord2 = (type2 == 'G' ? 0u : (type2 == 'P' ? 1u : 2u));
  value = 3u * ord1 + ord2;
}

aa_heur_pair::aa_heur_pair(const amino_acid &a1, const amino_acid &a2) {
  auto code1 = static_cast<aa_code>(a1);
  auto type1 =
      code1 == aa_code::GLY ? 'G' : (code1 == aa_code::PRO ? 'P' : 'X');

  auto code2 = static_cast<aa_code>(a2);
  auto type2 =
      code2 == aa_code::GLY ? 'G' : (code2 == aa_code::PRO ? 'P' : 'X');

  *this = aa_heur_pair(type1, type2);
}
aa_heur_pair::operator uint8_t() const { return value; }
