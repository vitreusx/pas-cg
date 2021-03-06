#include <cg/input/parameters.h>

namespace cg::input {

void parameters::load(ioxx::xyaml::node const &p) {
  if (auto pdb_node = p["pdb file"]; pdb_node && pdb_node["source"]) {
    pdb_source pdb_source;
    pdb_node >> pdb_source.file;

    if (auto deriv_node = pdb_node["contact deriv"]; deriv_node) {
      auto deriv_type = deriv_node.as<std::string>();
      if (deriv_type == "only from residues")
        pdb_source.deriv = pdb_file::contact_deriv::FROM_RESIDUES;
      else if (deriv_type == "from all atoms")
        pdb_source.deriv = pdb_file::contact_deriv::FROM_ATOMS;
    }

    if (auto ignore_node = pdb_node["ignore CRYST1"]; ignore_node)
      pdb_source.ignore_cryst1 = ignore_node.as<bool>();

    source = pdb_source;
  } else if (auto sf_node = p["seq file"]; sf_node && sf_node["source"]) {
    seq_file sf;
    sf_node["source"] >> sf;
    source = sf;
  }

  if (p["morph into SAW"])
    p["morph into SAW"] >> morph_into_saw;

  p["normalize mass"] >> normalize_mass;
  p["load native structure"] >> load_structure;
}
} // namespace cg::input