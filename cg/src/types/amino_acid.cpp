#include "types/amino_acid.h"
#include "utils/quantity.h"
#include <numeric>
#include <regex>
using namespace cg;

size_t std::hash<amino_acid>::operator()(const cg::amino_acid &aa) const {
  return std::hash<char>()(static_cast<char>((aa_code)aa));
}

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

void contact_limits::connect(ioxx::xyaml_proxy &proxy) {
  proxy["back"] & back;
  proxy["side-all"] & side_all;
  proxy["side-hydrophobic"] & side_hydrophobic;
  proxy["side-polar"] & side_polar;
}

void amino_acid_data::connect(ioxx::xyaml_proxy &proxy) {
  using namespace ioxx;

  if (proxy.loading()) {
    xyaml_embedded embed;
    proxy &embed;
    auto real_proxy = proxy(embed.node);

    std::unordered_map<std::string, double> def_atom_radii;
    for (auto const &entry : real_proxy["default atom radii"]) {
      auto name = entry.first.as<std::string>();
      auto radius = quantity(entry.second.as<std::string>()).in("A");
      def_atom_radii[name] = radius;
    }

    for (auto const &entry : real_proxy["amino acids"]) {
      auto name = entry.first.as<std::string>();
      auto data_proxy = real_proxy(entry.second);

      aa_data &cur_data = data[amino_acid(name)];
      cur_data.mass = quantity(data_proxy["mass"].as<std::string>()).in("amu");
      cur_data.radius =
          quantity(data_proxy["radius"].as<std::string>()).in("angstrom");

      auto atom_radii = def_atom_radii;
      if (auto alt = data_proxy["alt atom radii"]; alt) {
        for (auto const &alt_entry : alt) {
          auto atom_name = alt_entry.first.as<std::string>();
          auto alt_radius =
              quantity(alt_entry.second.as<std::string>()).in("A");
          atom_radii[atom_name] = alt_radius;
        }
      }

      for (auto const &back_atom : {"N", "CA", "C", "O", "OXT"}) {
        atom_data &atom_ = cur_data.atoms[back_atom];
        atom_.name = back_atom;
        atom_.radius = atom_radii[back_atom];
        atom_.backbone = true;
      }

      if (auto side_atoms_node = data_proxy["side"]; side_atoms_node) {
        for (auto const &side_atom_node : side_atoms_node) {
          auto side_atom = side_atom_node.as<std::string>();
          atom_data &atom_ = cur_data.atoms[side_atom];
          atom_.name = side_atom;
          atom_.radius = atom_radii[side_atom];
          atom_.backbone = false;
        }
      }

      cur_data.polarization = polarization_type::NONE;
      if (auto polarization_node = data_proxy["polarization"];
          polarization_node) {
        auto ptype = polarization_node.as<std::string>();
        if (ptype == "polar") {
          cur_data.polarization = polarization_type::POLAR;
        } else if (ptype == "hydrophobic") {
          cur_data.polarization = polarization_type::HYDROPHOBIC;
        }
      }

      cur_data.charge = 0.0;
      if (auto charge_node = data_proxy["charge"]; charge_node) {
        auto charge = quantity(charge_node.as<std::string>()).in("e");
        cur_data.charge = charge;
      }

      data_proxy["contact limits"] & cur_data.limits;
    }

    auto avg_mass = std::accumulate(
        data.begin(), data.end(), 0.0,
        [](auto const &sum, auto const &entry) -> auto {
          auto const &[name, res_data] = entry;
          return sum + res_data.mass;
        });
    avg_mass /= (double)data.size();

    for (auto &[name, res_data] : data) {
      res_data.mass = quantity(res_data.mass / avg_mass, "f77mass");
    }
  }
}

aa_data const &amino_acid_data::operator[](const amino_acid &aa) const {
  return data.at(aa);
}

bool cg::operator<(amino_acid const &aa1, amino_acid const &aa2) {
  return (uint8_t)aa1 < (uint8_t)aa2;
}

bool cg::operator<=(amino_acid const &aa1, amino_acid const &aa2) {
  return (uint8_t)aa1 <= (uint8_t)aa2;
}

bool cg::operator>(amino_acid const &aa1, amino_acid const &aa2) {
  return (uint8_t)aa1 > (uint8_t)aa2;
}

bool cg::operator>=(amino_acid const &aa1, amino_acid const &aa2) {
  return (uint8_t)aa1 >= (uint8_t)aa2;
}

bool cg::operator==(amino_acid const &aa1, amino_acid const &aa2) {
  return (uint8_t)aa1 == (uint8_t)aa2;
}

bool cg::operator!=(amino_acid const &aa1, amino_acid const &aa2) {
  return (uint8_t)aa1 != (uint8_t)aa2;
}