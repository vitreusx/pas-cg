#pragma once
#include "parameters.h"
#include <cg/input/pdb_file.h>
#include <cg/nat_cont/eval_forces.h>
#include <cg/pid/eval_forces.h>
#include <cg/qa/process_contacts.h>
#include <cg/simul/state.h>
#include <cg/types/amp.h>
#include <limits>

namespace cg::out {
class make_report {
public:
  std::filesystem::path prefix;
  real stats_every, struct_every;

public:
  report *rep;
  simul::state const *st;

  nat_cont::eval_forces const *nc;
  qa::process_contacts const *qa;
  pid::eval_forces const *pid;

public:
  void operator()() const;

private:
  real kinetic_energy() const;

  static real rmsd(vect::const_view<vec3r> const &orig_r,
                   vect::const_view<vec3r> const &cur_r,
                   vect::const_view<int> const &indices);

  static real gyration_radius(vect::const_view<vec3r> const &r,
                              vect::const_view<int> const &indices);

  static real asphericity(vect::const_view<vec3r> const &r,
                          vect::const_view<int> const &indices);

  void add_cur_scalars() const;
  void emit_out() const;

  void add_cur_snapshot() const;
  void emit_pdb() const;

  void add_map_data() const;
  void emit_map() const;
};
} // namespace cg::out