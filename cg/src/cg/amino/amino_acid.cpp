#include <cg/amino/amino_acid.h>
#include <cg/utils/quantity.h>
#include <regex>

namespace std {
size_t hash<cg::amino_acid>::operator()(const cg::amino_acid &aa) const {
  return std::hash<char>()(static_cast<char>((cg::aa_code)aa));
}
} // namespace std

namespace cg {

amino_acid::amino_acid(const aa_code &code) : code{code} {};

amino_acid::amino_acid(char letter) {
  static const std::unordered_map<char, aa_code> letter_to_code = {
      {'A', aa_code::ALA}, {'R', aa_code::ARG}, {'N', aa_code::ASN},
      {'D', aa_code::ASP}, {'C', aa_code::CYS}, {'E', aa_code::GLU},
      {'Q', aa_code::GLN}, {'G', aa_code::GLY}, {'H', aa_code::HIS},
      {'I', aa_code::ILE}, {'L', aa_code::LEU}, {'K', aa_code::LYS},
      {'M', aa_code::MET}, {'F', aa_code::PHE}, {'P', aa_code::PRO},
      {'S', aa_code::SER}, {'T', aa_code::THR}, {'W', aa_code::TRP},
      {'Y', aa_code::TYR}, {'V', aa_code::VAL}};

  code = letter_to_code.at(letter);
}

amino_acid::amino_acid(const std::string &name) {
  static const std::unordered_map<std::string, aa_code> name_to_code = {
      {"ALA", aa_code::ALA}, {"ARG", aa_code::ARG}, {"ASN", aa_code::ASN},
      {"ASP", aa_code::ASP}, {"CYS", aa_code::CYS}, {"GLU", aa_code::GLU},
      {"GLN", aa_code::GLN}, {"GLY", aa_code::GLY}, {"HIS", aa_code::HIS},
      {"ILE", aa_code::ILE}, {"LEU", aa_code::LEU}, {"LYS", aa_code::LYS},
      {"MET", aa_code::MET}, {"PHE", aa_code::PHE}, {"PRO", aa_code::PRO},
      {"SER", aa_code::SER}, {"THR", aa_code::THR}, {"TRP", aa_code::TRP},
      {"TYR", aa_code::TYR}, {"VAL", aa_code::VAL}};

  code = name_to_code.at(name);
}

char amino_acid::letter() const {
  static std::string letters = "ARNDCEQGHILKMFPSTWYV";
  return letters[static_cast<char>(code)];
}

std::string const &amino_acid::name() const {
  static std::vector<std::string> names = {
      "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
      "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};

  return names[static_cast<char>(code)];
}

std::vector<amino_acid> const &amino_acid::all() {
  struct amino_acids {
    std::vector<amino_acid> data;

    amino_acids() {
      data = std::vector<amino_acid>(amino_acid::NUM_TYPES);
      for (int i = 0; i < amino_acid::NUM_TYPES; ++i)
        data[i] = amino_acid(static_cast<aa_code>(i));
    }
  };

  static amino_acids amino_acids_;
  return amino_acids_.data;
}

std::ostream &operator<<(std::ostream &os, amino_acid const &aa) {
  os << aa.letter();
  return os;
}

std::istream &operator>>(std::istream &is, amino_acid &aa) {
  char letter;
  is >> letter;
  aa = amino_acid(letter);
  return is;
}

std::ostream &operator<<(std::ostream &os,
                         std::vector<amino_acid> const &acids) {
  std::string acids_str;
  for (size_t i = 0; i < acids.size(); ++i)
    acids_str[i] = acids[i].letter();
  os << acids_str;
  return os;
}

std::istream &operator>>(std::istream &is, std::vector<amino_acid> &acids) {
  std::string acids_str;
  is >> acids_str;
  acids_str = std::regex_replace(acids_str, std::regex("\\s+"), "");
  acids.resize(acids_str.size());
  for (size_t i = 0; i < acids.size(); ++i)
    acids[i] = amino_acid(acids_str[i]);
  return is;
}

bool operator<(amino_acid const &aa1, amino_acid const &aa2) {
  return (uint8_t)aa1 < (uint8_t)aa2;
}

bool operator<=(amino_acid const &aa1, amino_acid const &aa2) {
  return (uint8_t)aa1 <= (uint8_t)aa2;
}

bool operator>(amino_acid const &aa1, amino_acid const &aa2) {
  return (uint8_t)aa1 > (uint8_t)aa2;
}

bool operator>=(amino_acid const &aa1, amino_acid const &aa2) {
  return (uint8_t)aa1 >= (uint8_t)aa2;
}

bool operator==(amino_acid const &aa1, amino_acid const &aa2) {
  return (uint8_t)aa1 == (uint8_t)aa2;
}

bool operator!=(amino_acid const &aa1, amino_acid const &aa2) {
  return (uint8_t)aa1 != (uint8_t)aa2;
}
} // namespace cg