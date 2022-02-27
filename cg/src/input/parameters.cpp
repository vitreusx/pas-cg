#include "input/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::input;

void parameters::morph_into_saw_t::link(ioxx::xyaml::proxy &p) {
  p["perform"] & perform;
  p["bond distance"] & bond_distance;
  p["residue density"] & residue_density;
  p["infer simulation box"] & infer_box;
}

void parameters::load(ioxx::xyaml::node const &p) {
  auto type = p["type"].as<std::string>();
  if (type == "pdb file") {
    pdb_file file;
    p["source"] >> file;
    source = std::move(file);
  } else {
    seq_file file;
    p["source"] >> file;
    source = std::move(file);
  }

  if (p["morph into SAW"])
    p["morph into SAW"] >> morph_into_saw;
}