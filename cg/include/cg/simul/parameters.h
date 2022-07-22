#pragma once
#include <cg/afm/parameters.h>
#include <cg/amino/aa_data.h>
#include <cg/angles/parameters.h>
#include <cg/chir/parameters.h>
#include <cg/ckpt/parameters.h>
#include <cg/dh/parameters.h>
#include <cg/general/parameters.h>
#include <cg/input/parameters.h>
#include <cg/langevin/parameters.h>
#include <cg/local_rep/parameters.h>
#include <cg/nat_cont/parameters.h>
#include <cg/nl/parameters.h>
#include <cg/output/parameters.h>
#include <cg/pauli/parameters.h>
#include <cg/pbar/parameters.h>
#include <cg/pid/parameters.h>
#include <cg/qa/parameters.h>
#include <cg/sbox/parameters.h>
#include <cg/tether/parameters.h>

namespace cg::simul {
struct parameters {
  std::string repr;

  chir::parameters chir;
  dh::parameters dh;
  angles::parameters angles;
  lang::parameters lang;
  input::parameters input;
  nat_cont::parameters nat_cont;
  nl::parameters nl;
  pbar::parameters pbar;
  pid::parameters pid;
  tether::parameters tether;
  gen::parameters gen;
  amino_acid_data aa_data;
  qa::parameters qa;
  pauli::parameters pauli;
  out::parameters out;
  afm::parameters afm;
  ckpt::parameters ckpt;
  local_rep::parameters lrep;
  sbox::parameters sbox;

  void load(ioxx::xyaml::node const &p);
};

ioxx::xyaml::node defaults_yml();
} // namespace cg::simul