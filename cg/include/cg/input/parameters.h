#pragma once
#include "morph_into_saw.h"
#include "pdb_file.h"
#include "seq_file.h"
#include <cg/files/files.h>
#include <cg/utils/quantity.h>
#include <variant>

namespace cg::input {
struct parameters {
  struct pdb_source {
    pdb_file file;
    std::optional<pdb_file::contact_deriv> deriv;
    std::unordered_map<std::string, amino_acid> res_aliases;
    bool ignore_cryst1;
  };

  std::variant<seq_file, pdb_source> source;
  std::optional<morph_into_saw_t> morph_into_saw;
  bool normalize_mass, load_structure;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::input