#pragma once
#include "pdb_file.h"
#include "seq_file.h"
#include <cg/utils/quantity.h>
#include <ioxx/ioxx.h>
#include <variant>

namespace cg::input {
struct parameters {
  std::variant<seq_file, pdb_file> source;

  struct morph_into_saw_t {
    bool perform;
    std::optional<quantity> bond_distance;
    quantity residue_density;
    bool infer_box;
    void link(ioxx::xyaml::proxy &proxy);
  };
  std::optional<morph_into_saw_t> morph_into_saw;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::input