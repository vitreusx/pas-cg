#include "features/parameters.h"
using namespace cg;

parameters::parameters(std::filesystem::path const &path) {
  using namespace ioxx;
  auto paramfile = xyaml_node::from_path(path);
  auto proxy = xyaml_proxy(paramfile, node_proxy_mode::LOAD);
  proxy &*this;
}

void parameters::connect(ioxx::xyaml_proxy &p) {
  p["chirality"] & chir;
  p["complex native dihedrals"] & cnd;
  p["constant Debye-Hueckel"] & const_dh;
  p["force AFM"] & fafm;
  p["heurestic angles"] & heur_ang;
  p["heurestic dihedrals"] & heur_dih;
  p["input"] & input;
  p["native angles"] & nat_ang;
  p["native contacts"] & nat_cont;
  p["pseudo-improper dihedral"] & pid;
  p["quasi-adiabatic"] & qa;
  p["relative Debye-Hueckel"] & rel_dh;
  p["simple native dihedrals"] & snd;
  p["tether forces"] & tether;
  p["velocity AFM"] & vafm;
  p["Pauli exclusion"] & pauli;
  p["neighbor list"] & nl;
}