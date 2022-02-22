#pragma once
#include "chir/parameters.h"
#include "dh/const_dh/parameters.h"
#include "dh/rel_dh/parameters.h"
#include "force_afm/parameters.h"
#include "heur/heur_ang/parameters.h"
#include "heur/heur_dih/parameters.h"
#include "input/parameters.h"
#include "nat_ang/parameters.h"
#include "nat_cont/parameters.h"
#include "nat_dih/comp_nat_dih/parameters.h"
#include "nat_dih/simp_nat_dih/parameters.h"
#include "nl/parameters.h"
#include "pauli/parameters.h"
#include "pid/parameters.h"
#include "qa/parameters.h"
#include "tether/parameters.h"
#include "vel_afm/parameters.h"
#include <filesystem>

namespace cg {
struct parameters {
  parameters() = default;
  explicit parameters(std::filesystem::path const &path);

  chir::parameters chir;
  cnd::parameters cnd;
  const_dh::parameters const_dh;
  fafm::parameters fafm;
  heur_ang::parameters heur_ang;
  heur_dih::parameters heur_dih;
  input::parameters input;
  nat_ang::parameters nat_ang;
  nat_cont::parameters nat_cont;
  pid::parameters pid;
  qa::parameters qa;
  rel_dh::parameters rel_dh;
  snd::parameters snd;
  tether::parameters tether;
  vafm::parameters vafm;
  pauli::parameters pauli;
  nl::parameters nl;

  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg