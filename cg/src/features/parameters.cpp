#include "features/parameters.h"
using namespace cg;

void parameters::connect(ioxx::xyaml_proxy &p) {
  p["chirality"] & chir;
  p["complex native dihedral"] & cnd;
  p["constant Debye-Hueckel"] & const_dh;
  p["force AFM"] & fafm;
  p["heurestic angle"] & heur_ang;
  p["heurestic dihedral"] & heur_dih;
  p["input"] & input;
  p["native angle"] & nat_ang;
  p["native contacts"] & nat_cont;
  p["pseudo-improper dihedral"] & pid;
  p["quasi-adiabatic"] & qa;
  p["relative Debye-Hueckel"] & rel_dh;
  p["simple native dihedral"] & snd;
  p["tether"] & tether;
  p["velocity AFM"] & vafm;
}