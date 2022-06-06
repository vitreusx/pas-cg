#pragma once
#include "data.h"
#include "dynamics.h"
#include "parameters.h"
#include <fstream>
#include <iostream>
#include <mutex>
#include <omp.h>
#include <thread>

namespace cg::simul {
class state;

std::istream &operator>>(std::istream &is, state &st);
std::ostream &operator<<(std::ostream &os, state const &st);

class state {
public:
  explicit state(parameters &params);
  void verify_equal(state const &other) const;

  parameters params;

  bool did_simul_setup = false;
  bool is_running;
  real total_time, equil_time;
  rand_gen gen;
  void simul_setup();

  input::model orig_model;
  int num_res;
  void load_model();

  bool did_traj_setup = false;
  int traj_idx;
  void traj_setup();
  void finish_trajectory();

  input::model model;
  void morph_model();

  vect::vector<vec3r> orig_r, r;
  vect::vector<amino_acid> atype;
  compiled_aa_data comp_aa_data;
  cg::box<real> box;
  vect::vector<int> prev, next, chain_idx, seq_idx, chain_first, chain_last;
  input::model::res_map_t res_map;
  void compile_model();

  real t;
  dynamics dyn;
  void setup_dyn();

  out::report rep;
  real ckpt_last_t;
  void setup_output();

  vect::vector<real> mass_inv, mass_rsqrt;
  vect::vector<vec3r> v;
  vect::vector<vec3sr> y0, y1, y2, y3, y4, y5;
  solver_real true_t;
  real temperature;
  void setup_langevin();

  bool pbar_first_time;
  pbar::render::time_point_t pbar_start_clock, pbar_last_clock;
  void setup_pbar();

  nl::data nl;
  bool nl_required, nl_invalid, verify_first_time;
  real fixed_cutoff, max_cutoff, *max_cutoff_ptr, total_disp;
  vect::vector<int> res_cell_idx, reordered_idx, num_res_in_cell, cell_offset;
  vect::vector<nl::pair> all_pairs;
  void setup_nl();

  vect::vector<local_rep::pair> local_rep_pairs;
  void setup_local_rep();

  vect::vector<chir::chiral_quad> chir_quads;
  void setup_chir();

  vect::vector<tether::pair> tether_pairs;
  void setup_tether();

  vect::vector<nat_ang::nat_ang> native_angles;
  void setup_nat_ang();

  vect::vector<heur_ang::heur_ang> heur_angles;
  void setup_heur_ang();

  vect::vector<nat_dih> native_dihedrals;
  void setup_nat_dih();

  vect::vector<heur_dih::heur_dih> heur_dihedrals;
  void setup_heur_dih();

  vect::vector<pauli::pair> pauli_pairs;
  real pauli_cutoff;
  void setup_pauli();

  vect::vector<nat_cont::nat_cont> all_native_contacts, cur_native_contacts;
  vect::vector<nl::exclusion> nat_cont_excl;
  real nat_cont_cutoff;
  void setup_nat_cont();

  vect::vector<dh::pair> dh_pairs;
  void setup_dh();

  vect::set<qa::free_pair> qa_free_pairs;
  vect::vector<qa::candidate> qa_candidates;
  vect::set<qa::contact> qa_contacts;
  vect::vector<sync_data> sync_values;
  vect::vector<vec3r> n, h;
  vect::vector<qa::cys_neigh> qa_cys_neigh;
  bool ss_spec_crit;
  vect::vector<int> neigh_count, cys_indices, qa_removed;
  vect::vector<bool> part_of_ssbond;
  int num_qa_contacts;
  void setup_qa();

  vect::vector<pid::bundle> pid_bundles;
  vect::vector<sink_lj> ss_ljs;
  void setup_pid();

public:
  bool post_equil = false;
  bool did_post_equil_setup = false;
  void post_equil_setup();

  afm::compiled_tips afm_tips;
  void setup_afm();

public:
  friend std::istream &operator>>(std::istream &is, state &st);
  friend std::ostream &operator<<(std::ostream &os, state const &st);
};
} // namespace cg::simul