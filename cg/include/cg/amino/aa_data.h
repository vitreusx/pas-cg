#pragma once
#include <cg/amino/amino_acid.h>
#include <cg/files/files.h>
#include <cg/utils/quantity.h>

namespace cg {

struct atom_data {
  std::string name;
  quantity radius;
  bool backbone;
};

struct contact_limits {
  int back, side_all, side_hydrophobic, side_polar;
  void link(ioxx::xyaml::proxy &proxy);
};

struct aa_data {
  quantity mass, radius;
  std::unordered_map<std::string, atom_data> atoms;
  polarization_type polarization;
  quantity charge;
  contact_limits limits;

  atom_data const& for_atom(std::string const& atom_name) const;
};

class amino_acid_data {
public:
  amino_acid_data() = default;
  std::unordered_map<amino_acid, aa_data> data;
  aa_data const &operator[](amino_acid const &aa) const;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg
