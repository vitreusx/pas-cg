#pragma once
#include <cg/afm/force/eval_forces.h>
#include <cg/afm/vel/eval_forces.h>
#include <cg/chir/eval_forces.h>
#include <cg/dh/const/eval_forces.h>
#include <cg/dh/rel/eval_forces.h>
#include <cg/dh/update_pairs.h>
#include <cg/heur/ang/eval_forces.h>
#include <cg/heur/dih/eval_forces.h>
#include <cg/langevin/legacy_step.h>
#include <cg/langevin/step.h>
#include <cg/nat_ang/eval_forces.h>
#include <cg/nat_cont/eval_forces.h>
#include <cg/nat_cont/report.h>
#include <cg/nat_cont/update_contacts.h>
#include <cg/nat_dih/complex/eval_forces.h>
#include <cg/nat_dih/simple/eval_forces.h>
#include <cg/nl/cell_update.h>
#include <cg/nl/legacy_update.h>
#include <cg/nl/verify.h>
#include <cg/output/export_pdb.h>
#include <cg/output/gyration.h>
#include <cg/output/make_report.h>
#include <cg/output/rmsd.h>
#include <cg/output/stats.h>
#include <cg/output/structure.h>
#include <cg/pauli/eval_forces.h>
#include <cg/pauli/update_pairs.h>
#include <cg/pbar/render.h>
#include <cg/pid/eval_forces.h>
#include <cg/pid/update_bundles.h>
#include <cg/qa/count_cys_neigh.h>
#include <cg/qa/finish_processing.h>
#include <cg/qa/prepare_nh.h>
#include <cg/qa/process_contacts.h>
#include <cg/qa/report.h>
#include <cg/qa/sift_candidates.h>
#include <cg/qa/update_free_pairs.h>
#include <cg/tether/eval_forces.h>
#include <cg/utils/random.h>

namespace cg::simul {
class kernels {
public:
  chir::eval_forces eval_chir_forces;
  const_dh::eval_forces eval_const_dh_forces;
  rel_dh::eval_forces eval_rel_dh_forces;
  afm::force::eval_forces eval_force_afm_forces;
  heur_ang::eval_forces eval_heur_ang_forces;
  heur_dih::eval_forces eval_heur_dih_forces;
  nat_ang::eval_forces eval_nat_ang_forces;
  nat_cont::eval_forces eval_nat_cont_forces;
  cnd::eval_forces eval_cnd_forces;
  snd::eval_forces eval_snd_forces;
  pauli::eval_forces eval_pauli_forces;
  pid::eval_forces eval_pid_forces;
  tether::eval_forces eval_tether_forces;
  afm::vel::eval_forces eval_vel_afm_forces;

  lang::step lang_step;
  lang::legacy_step lang_legacy_step;
  nl::legacy_update nl_legacy;
  nl::cell_update nl_cell;
  nl::verify nl_verify;

  dh::update_pairs update_dh_pairs;
  nat_cont::update_contacts update_nat_contacts;
  pauli::update_pairs update_pauli_pairs;
  pid::update_bundles update_pid_bundles;
  qa::update_free_pairs update_qa_pairs;
  qa::update_cys_neigh update_cys_neigh;

  qa::prepare_nh prepare_nh;
  qa::count_cys_neigh count_cys_neigh;
  qa::sift_candidates sift_qa_candidates;
  qa::process_contacts process_qa_contacts;
  qa::finish_processing qa_finish_processing;

  pbar::render render_pbar;
  out::make_report make_report;
  out::export_pdb export_pdb;
  out::add_stats add_stats;
  out::add_structure add_structure;
  qa::report_qa_stuff report_qa_stuff;
  nat_cont::report_stuff report_nc_stuff;
  out::report_gyration_stuff report_gyr;
  out::compute_rmsd compute_rmsd;
};
} // namespace cg::simul