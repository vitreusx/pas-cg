#pragma once
#include <cg/afm/parameters.h>
#include <cg/amino/aa_data.h>
#include <cg/base_forces/lj_variants.h>
#include <cg/chir/parameters.h>
#include <cg/ckpt/parameters.h>
#include <cg/dh/const/parameters.h>
#include <cg/dh/rel/parameters.h>
#include <cg/general/parameters.h>
#include <cg/heur/ang/parameters.h>
#include <cg/heur/dih/parameters.h>
#include <cg/input/parameters.h>
#include <cg/langevin/parameters.h>
#include <cg/local_rep/parameters.h>
#include <cg/nat_ang/parameters.h>
#include <cg/nat_cont/parameters.h>
#include <cg/nat_dih/complex/parameters.h>
#include <cg/nat_dih/simple/parameters.h>
#include <cg/nl/parameters.h>
#include <cg/output/parameters.h>
#include <cg/pauli/parameters.h>
#include <cg/pbar/parameters.h>
#include <cg/pid/parameters.h>
#include <cg/qa/parameters.h>
#include <cg/tether/parameters.h>

namespace cg::simul {
struct parameters {
  std::string repr;

  chir::parameters chir;
  const_dh::parameters const_dh;
  rel_dh::parameters rel_dh;
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
  gen::parameters gen;
  amino_acid_data aa_data;
  lj_variants lj_vars;
  qa::parameters qa;
  pauli::parameters pauli;
  out::parameters out;
  afm::parameters afm;
  ckpt::parameters ckpt;
  local_rep::parameters lrep;

  void load(ioxx::xyaml::node const &p);
};

ioxx::xyaml::node defaults_yml();
} // namespace cg::simul