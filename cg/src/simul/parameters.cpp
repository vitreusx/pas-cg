#include "simul/parameters.h"
using namespace cg::simul;

void parameters::connect(ioxx::xyaml_proxy &proxy) {
  proxy["general"] & gen;
  proxy["amino acid data"] & aa_data;
  proxy["input"] & input;
  proxy["langevin"] & lang;
  proxy["neighbor list"] & nl;
  proxy["heurestic angles"] & heur_ang;
  proxy["native angles"] & nat_ang;
  proxy["complex native dihedrals"] & cnd;
  proxy["simple native dihedrals"] & snd;
  proxy["heurestic dihedrals"] & heur_dih;
  proxy["constant Debye-Hueckel"] & const_dh;
  proxy["relative Debye-Hueckel"] & rel_dh;
  proxy["lj force variants"] & lj_variants;
  proxy["pseudo-improper dihedral"] & pid;
  proxy["quasi-adiabatic"] & qa;
  proxy["chirality"] & chir;
  proxy["native contacts"] & nat_cont;
  proxy["Pauli exclusion"] & pauli;
  proxy["tether forces"] & tether;
  proxy["velocity AFM"] & vafm;
  proxy["force AFM"] & fafm;
  proxy["progress bar"] & pbar;
}