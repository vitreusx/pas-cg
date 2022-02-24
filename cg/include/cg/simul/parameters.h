#pragma once
#include <cg/amino/aa_data.h>
#include <cg/base_forces/lj_variants.h>
#include <cg/chir/parameters.h>
#include <cg/dh/const/parameters.h>
#include <cg/dh/rel/parameters.h>
#include <cg/force_afm/parameters.h>
#include <cg/general/parameters.h>
#include <cg/heur/ang/parameters.h>
#include <cg/heur/dih/parameters.h>
#include <cg/input/parameters.h>
#include <cg/langevin/parameters.h>
#include <cg/nat_ang/parameters.h>
#include <cg/nat_cont/parameters.h>
#include <cg/nat_dih/complex/parameters.h>
#include <cg/nat_dih/simple/parameters.h>
#include <cg/nl/parameters.h>
#include <cg/pauli/parameters.h>
#include <cg/pbar/parameters.h>
#include <cg/pid/parameters.h>
#include <cg/qa/parameters.h>
#include <cg/tether/parameters.h>
#include <cg/vel_afm/parameters.h>

namespace cg::simul {
struct parameters {
  chir::parameters chir;
  const_dh::parameters const_dh;
  rel_dh::parameters rel_dh;
  fafm::parameters fafm;
  heur_ang::parameters heur_ang;
  heur_dih::parameters heur_dih;
  lang::parameters lang;
  input::parameters input;
  nat_ang::parameters nat_ang;
  nat_cont::parameters nat_cont;
  cnd::parameters cnd;
  snd::parameters snd;
  nl::parameters nl;
  pbar::parameters pbar;
  pid::parameters pid;
  tether::parameters tether;
  vafm::parameters vafm;
  gen::parameters gen;
  amino_acid_data aa_data;
  lj_variants lj_variants;
  qa::parameters qa;
  pauli::parameters pauli;

  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg::simul