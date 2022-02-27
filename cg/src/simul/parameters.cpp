#include "simul/parameters.h"
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
  from["lj force variants"] >> lj_variants;
  from["pseudo-improper dihedral"] >> pid;
  from["quasi-adiabatic"] >> qa;
  from["chirality"] >> chir;
  from["native contacts"] >> nat_cont;
  from["Pauli exclusion"] >> pauli;
  from["tether forces"] >> tether;
  from["velocity AFM"] >> vafm;
  from["force AFM"] >> fafm;
  from["progress bar"] >> pbar;
}