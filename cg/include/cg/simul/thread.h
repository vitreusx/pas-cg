#pragma once
#include "state.h"

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <cfenv>

namespace cg::simul {
class thread {
public:
  explicit thread(state &st);
  thread(thread const &) = delete;

  void main();

  void adjust_scenario();
  void pre_eval();
  void eval_forces();
  void post_eval();
  void fix_nl();

public:
  state &st;
  parameters const &params;

public:
  bool did_overall_setup = false;
  void overall_setup();

  rand_gen gen;
  void setup_gen();

  dynamics dyn;
  void setup_dyn();

  std::vector<out::hook const *> hooks;
  out::make_report make_report;
  out::export_pdb export_pdb;
  out::add_stats add_stats;
  out::add_structure add_structure;
  qa::report_qa_stuff report_qa_stuff;
  nat_cont::report_stuff report_nc_stuff;
  out::report_gyration_stuff report_gyr;
  out::compute_rmsd compute_rmsd;
  void setup_output();

  lang::step lang_step;
  lang::legacy_step lang_legacy_step;
  void setup_langevin();

  pbar::render render_pbar;
  void setup_pbar();

  nl::legacy_update nl_legacy;
  nl::cell_update nl_cell;
  nl::verify nl_verify;
  void setup_nl();

  chir::eval_forces eval_chir_forces;
  void setup_chir();

  tether::eval_forces eval_tether_forces;
  void setup_tether();

  nat_ang::eval_forces eval_nat_ang_forces;
  void setup_nat_ang();

  heur_ang::eval_forces eval_heur_ang_forces;
  void setup_heur_ang();

  cnd::eval_forces eval_cnd_forces;
  snd::eval_forces eval_snd_forces;
  void setup_nat_dih();

  heur_dih::eval_forces eval_heur_dih_forces;
  void setup_heur_dih();

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
  qa::count_cys_neigh count_cys_neigh;
  qa::sift_candidates sift_qa_candidates;
  qa::process_contacts process_qa_contacts;
  qa::finish_processing qa_finish_processing;
  void setup_qa();

  pid::eval_forces eval_pid_forces;
  pid::update_bundles update_pid_bundles;
  void setup_pid();

public:
  bool did_post_equil_setup = false;
  void post_equil_setup();

  afm::force::eval_forces eval_force_afm_forces;
  afm::vel::eval_forces eval_vel_afm_forces;
  afm::report_stats report_afm_stats;
  void setup_afm();
};
} // namespace cg::simul