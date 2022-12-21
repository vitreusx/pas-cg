#pragma once
#include "kernels.h"
#include "state.h"

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <cfenv>

namespace cg::simul {
class thread;

class thread_team {
public:
  explicit thread_team(state &st);
  thread &fork();

public:
  state *st;
  int num_threads;
  std::vector<std::unique_ptr<thread>> threads;
  std::vector<vect::const_view<vec3r>> forces;
};

class thread {
public:
  explicit thread(thread_team &team);

  thread(thread const &) = delete;
  thread &operator=(thread const &) = delete;

  thread(thread &&) = delete;
  thread &operator=(thread &&) = delete;

  void main();

public:
  int tid;
  thread_team *team;

public:
  void step();
  void simul_init_step();
  void traj_init_step();
  void pull_release_step();
  void equil_step();
  void squeezing_step();
  void rest_after_squeezing_step();
  void find_force_min_step();
  void rest_after_force_min_step();
  void max_amplitude_step();
  void rest_after_max_amp_step();
  void oscillations_step();
  void rest_after_oscillations_step();
  void freeform_step();
  void traj_end_step();
  void simul_end_step();

  void init_kernels();

  void advance_by_step();
  void pre_eval_async();
  void fix_def_nl_async();
  void fix_lr_nl_async();
  void eval_forces();
  void post_eval_async();

public:
  int loop_idx;
  state *st;
  parameters const *params;

  set_of_divisibles eval_divs;

public:
  rand_gen gen;
  void setup_gen();

  dynamics dyn;
  void setup_dyn();

  out::make_report make_report;
  ckpt::make_checkpoint make_checkpoint;
  out::print_raw_data print_raw_data;
  void setup_output();

  lang::step lang_step;
  void setup_langevin();

  pbar::render render_pbar;
  void setup_pbar();

  struct nl_stuff {
    nl::legacy_update legacy;
    nl::cell_update cell;
    nl::verify verify;
  };

  nl_stuff def_nl, long_range_nl;
  void setup_nl();

  local_rep::eval_forces eval_lrep_forces;
  void setup_local_rep();

  chir::eval_forces eval_chir_forces;
  void setup_chir();

  tether::eval_forces eval_tether_forces;
  void setup_tether();

  nat_ang::eval_forces eval_nat_ang_forces;
  heur_ang::eval_forces eval_heur_ang_forces;
  cnd::eval_forces eval_cnd_forces;
  snd::eval_forces eval_snd_forces;
  heur_dih::eval_forces eval_heur_dih_forces;
  void setup_angles();

  pauli::update_pairs update_pauli_pairs;
  pauli::eval_forces eval_pauli_forces;
  void setup_pauli();

  nat_cont::update_contacts update_nat_contacts;
  nat_cont::eval_forces eval_nat_cont_forces;
  void setup_nat_cont();

  dh::update_pairs update_dh_pairs;
  const_dh::eval_forces eval_const_dh_forces;
  rel_dh::eval_forces eval_rel_dh_forces;
  void setup_dh();

  qa::update_free_pairs update_qa_pairs;
  qa::update_cys_neigh update_cys_neigh;
  qa::prepare_nh prepare_nh;
  qa::reset_cys_neigh reset_cys_neigh;
  qa::count_cys_neigh count_cys_neigh;
  qa::loop_over_candidates qa_loop_over_candidates;
  qa::process_contacts process_qa_contacts;
  qa::finish_processing qa_finish_processing;
  void setup_qa();

  pid::eval_forces eval_pid_forces;
  pid::eval_forces::fast_version_t eval_pid_fast;
  pid::update_bundles update_pid_bundles;
  void setup_pid();

  afm::vel::eval_forces eval_vel_afm_forces;
  afm::force::eval_forces eval_force_afm_forces;
  void setup_afm();

  wall::solid::eval_forces eval_solid_wall_forces;
  void setup_solid_walls();

  wall::harmonic::eval_free hw_eval_free;
  wall::harmonic::eval_connected hw_eval_conn;
  void setup_harmonic_walls();

  wall::lj::sift_free ljw_sift_free;
  wall::lj::eval_connected ljw_eval_conn;
  wall::lj::process_candidates ljw_proc_cand;
  void setup_lj_walls();

  bool log_wall_forces_enabled;
  wall::log_forces log_wall_forces;
  void setup_walls();
};
} // namespace cg::simul