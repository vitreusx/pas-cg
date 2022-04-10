#include <cg/simul/default.yml.h>
#include <cg/simul/parameters.h>
namespace cg::simul {

void parameters::link(ioxx::xyaml::proxy &p) {
  p["general"] & gen;
  p["amino acid data"] & aa_data;
  p["input"] & input;
  p["langevin"] & lang;
  p["neighbor list"] & nl;
  p["heurestic angles"] & heur_ang;
  p["native angles"] & nat_ang;
  p["complex native dihedrals"] & cnd;
  p["simple native dihedrals"] & snd;
  p["heurestic dihedrals"] & heur_dih;
  p["constant Debye-Hueckel"] & const_dh;
  p["relative Debye-Hueckel"] & rel_dh;
  p["lj force variants"] & lj_vars;
  p["pseudo-improper dihedral"] & pid;
  p["quasi-adiabatic"] & qa;
  p["chirality"] & chir;
  p["native contacts"] & nat_cont;
  p["Pauli exclusion"] & pauli;
  p["tether forces"] & tether;
  p["AFM"] & afm;
  p["progress bar"] & pbar;
  p["output"] & out;

  if (p.is_loading) {
    if (gen.debug_mode.determinism) {
      out.enabled = false;
      pbar.enabled = false;
    }
  }
}

ioxx::xyaml::node defaults_yml() {
  auto default_yml_str =
      std::string((const char *)default_yml, default_yml_len);
  auto default_yml_node = YAML::Load(default_yml_str.c_str());
  return ioxx::xyaml::node(default_yml_node);
}
} // namespace cg::simul