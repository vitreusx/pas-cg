#include "features/input/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::input;

void parameters::morph_into_saw_t::connect(ioxx::xyaml_proxy &p) {
  p["bond distance"] & bond_distance;
  p["residue density"] & residue_density;
  p["infer box"] & infer_box;
}

void parameters::connect(ioxx::xyaml_proxy &p) {
  if (p["seq file"]) {
    seq_file sf;
    p["seq file"] & sf.path;
    source = sf;
  } else if (p["pdb file"]) {
    pdb_file pf;
    p["pdb file"] & pf.path;
    source = pf;
  }

  if (p["morph into SAW"])
    p["morph into SAW"] & morph_into_saw;
}