#pragma once
#include "data.h"
#include "dynamics.h"
#include "parameters.h"
#include <cg/types/avg.h>
#include <fstream>
#include <iostream>
#include <mutex>
#include <omp.h>
#include <thread>

namespace cg::simul {
enum class phase {
  SIMUL_INIT,
  TRAJ_INIT,
  PULL_RELEASE,
  EQUIL,
  SQUEEZING,
  REST_AFTER_SQUEEZING,
  FIND_FORCE_MIN,
  REST_AFTER_FORCE_MIN,
  MAX_AMPLITUDE,
  REST_AFTER_MAX_AMP,
  OSCILLATIONS,
  REST_AFTER_OSCILLATIONS,
  FREEFORM,
  TRAJ_END,
  SIMUL_END,
  PROG_END
};

enum side {
  POS_X = 0,
  NEG_X = 1,
  POS_Y = 2,
  NEG_Y = 3,
  POS_Z = 4,
  NEG_Z = 5
};

enum axis {
  X = 0,
  Y = 1,
  Z = 2
};

class state {
public:
  explicit state(parameters &params);
  void verify_equal(state const &other) const;

  parameters params;

  phase cur_phase = phase::SIMUL_INIT;
  int traj_idx, step_idx;
  real until, since, amplitude, displacement, omega;
  bool shear;

  rand_gen gen;
  void simul_setup();

  input::model orig_model;
  int num_res;
  real bond;
  void load_model();

  void traj_init();

  input::model model;
  void morph_model();

  vect::vector<vec3r> orig_r, r;
  vect::vector<amino_acid> atype;
  compiled_aa_data comp_aa_data;
  sbox::pbc<real> pbc;
  sbox::box<real> box;
  vect::vector<int> prev, next, chain_idx, seq_idx, chain_first, chain_last;
  input::model::res_map_t res_map;
  void compile_model();

  real t;
  dynamics dyn;
  void setup_dyn();

  bool out_enabled, ckpt_enabled, dump_enabled;
  out::report rep;
  real ckpt_last_t;
  void simul_setup_output();
  void setup_output();

  bool lang_enabled;
  vect::vector<real> mass_inv, mass_rsqrt;
  vect::vector<vec3r> v, noise;
  vect::vector<vec3sr> y0, y1, y2, y3, y4, y5;
  solver_real true_t;
  real temperature;
  void setup_langevin();

  bool pbar_enabled, pbar_first_time;
  pbar::render::time_point_t pbar_start_clock, pbar_last_clock;
  void setup_pbar();

  nl::data nl;
  bool nl_required, nl_invalid, verify_first_time;
  real total_disp;
  vect::vector<nl::pair> all_pairs;
  void setup_nl();

  bool lrep_enabled;
  vect::vector<local_rep::pair> local_rep_pairs;
  void setup_local_rep();

  bool chir_enabled;
  vect::vector<chir::chiral_quad> chir_quads;
  void setup_chir();

  bool tether_enabled;
  vect::vector<tether::pair> tether_pairs;
  void setup_tether();

  bool nat_ang_enabled;
  vect::vector<nat_ang::nat_ang> native_angles;
  void setup_nat_ang();

  bool heur_ang_enabled;
  vect::vector<heur_ang::heur_ang> heur_angles;
  void setup_heur_ang();

  bool nat_dih_enabled;
  vect::vector<nat_dih> native_dihedrals;
  void setup_nat_dih();

  bool heur_dih_enabled;
  vect::vector<heur_dih::heur_dih> heur_dihedrals;
  void setup_heur_dih();

  bool pauli_enabled;
  vect::vector<pauli::pair> pauli_pairs;
  void setup_pauli();

  bool nat_cont_enabled;
  vect::vector<nat_cont::nat_cont> all_native_contacts, cur_native_contacts;
  vect::vector<nl::exclusion> nat_cont_excl;
  int num_changed;
  vect::vector<vect::vector<real>> nc_times;
  vect::vector<real> nc_unfold_times;
  void setup_nat_cont();

  bool dh_enabled;
  vect::vector<dh::pair> dh_pairs;
  void setup_dh();

  bool qa_enabled;
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

  bool pid_enabled;
  vect::vector<pid::bundle> pid_bundles;
  void setup_pid();

  bool harmonic_walls_enabled;
  vect::vector<wall::harmonic::connection> harmonic_conns;
  vect::vector<wall::harmonic::wall> harmonic_walls;

  bool solid_walls_enabled;
  vect::vector<wall::solid::wall> solid_walls;

  bool lj_walls_enabled;
  vect::vector<int> ljw_removed;
  vect::set<wall::lj::connection> ljw_conns;
  vect::vector<wall::lj::candidate> ljw_candidates;
  vect::vector<wall::lj::wall> lj_walls;

  vect::vector<bool> is_connected_to_wall;
  std::string wall_type[3];
  bool pbc_on[3];
  vect::vector<wall::gen_wall *> walls;
  moving_avg<real, real> avg_z_force;

  void setup_walls();
  void adjust_wall_pos(vec3r size_change, vec3r translation);
  void reinit_wall_values();

  bool trajectory_should_end() const;

  bool afm_enabled;
  vect::vector<afm::vel::tip> vel_afm_tips;
  vect::vector<afm::force::tip> force_afm_tips;
  void setup_afm();

public:
  friend std::istream &operator>>(std::istream &is, state &st);
  friend std::ostream &operator<<(std::ostream &os, state const &st);
};
} // namespace cg::simul