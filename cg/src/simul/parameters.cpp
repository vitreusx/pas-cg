#include "simul/parameters.h"
#include "default.yml.h"
using namespace cg::simul;

void parameters::load(ioxx::xyaml::node const &from) {
  from["general"] >> gen;
  from["amino acid data"] >> aa_data;
  from["input"] >> input;
  from["langevin"] >> lang;
  from["neighbor list"] >> nl;
  from["heurestic angles"] >> heur_ang;
  from["native angles"] >> nat_ang;
  from["complex native dihedrals"] >> cnd;
  from["simple native dihedrals"] >> snd;
  from["heurestic dihedrals"] >> heur_dih;
  from["constant Debye-Hueckel"] >> const_dh;
  from["relative Debye-Hueckel"] >> rel_dh;
  from["lj force variants"] >> lj_vars;
  from["pseudo-improper dihedral"] >> pid;
  from["quasi-adiabatic"] >> qa;
  from["chirality"] >> chir;
  from["native contacts"] >> nat_cont;
  from["Pauli exclusion"] >> pauli;
  from["tether forces"] >> tether;
  from["AFM"] >> afm;
  from["progress bar"] >> pbar;
  from["output"] >> out;
}

ioxx::xyaml::node cg::simul::defaults_yml() {
  auto default_yml_str =
      std::string((const char *)default_yml, default_yml_len);
  auto default_yml_node = YAML::Load(default_yml_str.c_str());
  return ioxx::xyaml::node(default_yml_node);
}