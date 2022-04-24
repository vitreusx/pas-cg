#pragma once
#include "parameters.h"
#include <cg/input/pdb_file.h>
#include <cg/nat_cont/eval_forces.h>
#include <cg/pid/eval_forces.h>
#include <cg/qa/process_contacts.h>
#include <cg/types/amp.h>
#include <limits>

namespace cg::out {
class report {
public:
  real stats_last, struct_last;
  int model_serial;
  pdb_file full_pdb;

  struct traj_data {
    int traj_idx;
    ioxx::table::table scalars;
  };
  std::map<int, traj_data> traj;

  void simul_init();
  void traj_init(int traj_idx);
};

struct gyration {
  explicit gyration(nitro::const_view<vec3r> const &r,
                    nitro::const_view<int> const &indices);
  vec3r lambda;

  real radius() const;
  real asphericity() const;
};

class make_report {
public:
  std::filesystem::path prefix;
  real stats_every, struct_every;

public:
  report *rep;

  real const *t, *V;
  input::model const *model;
  input::model::res_map_t const *res_map;
  nitro::const_view<vec3r> r, orig_r, v;
  nitro::const_view<real> mass;
  nitro::const_view<amino_acid> atype;
  nitro::const_view<int> chain_first, chain_last;
  int *traj_idx;

  nat_cont::eval_forces const *nc;
  qa::process_contacts const *qa;
  pid::eval_forces const *pid;

public:
  void operator()() const;

private:
  real kinetic_energy() const;
  real rmsd(nitro::const_view<int> const &indices) const;

  void add_cur_scalars() const;
  void emit_out() const;

  void add_cur_model() const;
  void emit_pdb() const;

  void emit_map() const;
};
} // namespace cg::out