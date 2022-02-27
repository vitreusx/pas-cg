#pragma once
#include <string>
#include <unordered_map>
#include <vector>

namespace cg {
class amino_acid;
}

std::ostream &operator<<(std::ostream &os, cg::amino_acid const &aa);
std::istream &operator>>(std::istream &is, cg::amino_acid &aa);
std::ostream &operator<<(std::ostream &os,
                         std::vector<cg::amino_acid> const &acids);
std::istream &operator>>(std::istream &is, std::vector<cg::amino_acid> &acids);

namespace cg {
enum aa_code : char {
  ALA,
  ARG,
  ASN,
  ASP,
  CYS,
  GLU,
  GLN,
  GLY,
  HIS,
  ILE,
  LEU,
  LYS,
  MET,
  PHE,
  PRO,
  SER,
  THR,
  TRP,
  TYR,
  VAL
};

enum polarization_type : char { NONE, POLAR, HYDROPHOBIC };

class amino_acid {
public:
  amino_acid() = default;
  amino_acid(aa_code const &code);
  explicit amino_acid(char letter);
  explicit amino_acid(std::string const &name);

  inline operator aa_code() const { return code; };
  inline operator uint8_t() const { return static_cast<uint8_t>(code); }
  char letter() const;
  std::string const &name() const;

  static constexpr int NUM_TYPES = 20;
  static std::vector<amino_acid> const &all();

  friend std::ostream & ::operator<<(std::ostream &os, amino_acid const &aa);
  friend std::istream & ::operator>>(std::istream &is, amino_acid &aa);
  friend std::ostream & ::operator<<(std::ostream &os,
                                     std::vector<amino_acid> const &acids);
  friend std::istream & ::operator>>(std::istream &is,
                                     std::vector<amino_acid> &acids);

private:
  aa_code code;
};

bool operator<(amino_acid const &aa1, amino_acid const &aa2);
bool operator<=(amino_acid const &aa1, amino_acid const &aa2);
bool operator>(amino_acid const &aa1, amino_acid const &aa2);
bool operator>=(amino_acid const &aa1, amino_acid const &aa2);
bool operator==(amino_acid const &aa1, amino_acid const &aa2);
bool operator!=(amino_acid const &aa1, amino_acid const &aa2);

} // namespace cg

namespace std {
template <> struct hash<cg::amino_acid> {
  size_t operator()(cg::amino_acid const &aa) const;
};
} // namespace std
