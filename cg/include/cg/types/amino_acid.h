#pragma once
#include <ioxx/xyaml.h>
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

} // namespace cg

namespace std {
template <> struct hash<cg::amino_acid> {
  size_t operator()(cg::amino_acid const &aa) const;
};
} // namespace std

namespace cg {

struct atom_data {
  std::string name;
  double radius;
  bool backbone;
};

struct contact_limits {
  int back, side_all, side_hydrophobic, side_polar;
  void connect(ioxx::xyaml_proxy &proxy);
};

struct aa_data {
  double mass, radius;
  std::unordered_map<std::string, atom_data> atoms;
  polarization_type polarization;
  double charge;
  contact_limits limits;
};

class amino_acid_data {
public:
  amino_acid_data() = default;
  void connect(ioxx::xyaml_proxy &proxy);

  std::unordered_map<amino_acid, aa_data> data;
  aa_data const &operator[](amino_acid const &aa) const;
};
} // namespace cg
