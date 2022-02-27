#pragma once
#include "thread_state.h"
#include <cg/chir/eval_forces.h>
#include <cg/dh/const/eval_forces.h>
#include <cg/dh/rel/eval_forces.h>
#include <cg/dh/update_pairs.h>
#include <cg/force_afm/eval_forces.h>
#include <cg/heur/ang/eval_forces.h>
#include <cg/heur/dih/eval_forces.h>
#include <cg/langevin/step.h>
#include <cg/nat_ang/eval_forces.h>
#include <cg/nat_cont/eval_forces.h>
#include <cg/nat_cont/update_contacts.h>
#include <cg/nat_dih/complex/eval_forces.h>
#include <cg/nat_dih/simple/eval_forces.h>
#include <cg/nl/legacy_update.h>
#include <cg/nl/verify.h>
#include <cg/pauli/eval_forces.h>
#include <cg/pauli/update_pairs.h>
#include <cg/pbar/render.h>
#include <cg/pid/eval_forces.h>
#include <cg/pid/update_bundles.h>
#include <cg/qa/prepare_nh.h>
#include <cg/qa/process_candidates.h>
#include <cg/qa/process_contacts.h>
#include <cg/qa/sift_candidates.h>
#include <cg/qa/update_free_pairs.h>
#include <cg/tether/eval_forces.h>
#include <cg/utils/random.h>
#include <cg/vel_afm/eval_forces.h>
#include <cg/nl/cell_update.h>

namespace cg::simul {
class kernels {
public:
  chir::eval_forces eval_chir_forces;
  const_dh::eval_forces eval_const_dh_forces;
  rel_dh::eval_forces eval_rel_dh_forces;
  fafm::eval_forces eval_fafm_forces;
  heur_ang::eval_forces eval_heur_ang_forces;
  heur_dih::eval_forces eval_heur_dih_forces;
  nat_ang::eval_forces eval_nat_ang_forces;
  nat_cont::eval_forces eval_nat_cont_forces;
  cnd::eval_forces eval_cnd_forces;
  snd::eval_forces eval_snd_forces;
  pauli::eval_forces eval_pauli_forces;
  pid::eval_forces eval_pid_forces;
  qa::process_contacts process_qa_contacts;
  qa::sift_candidates sift_qa_candidates;
  tether::eval_forces eval_tether_forces;
  vafm::eval_forces eval_vafm_forces;

  lang::step lang_step;
  nl::legacy_update nl_legacy;
  nl::cell_update nl_cell;
  nl::verify nl_verify;

  dh::update_pairs update_dh_pairs;
  nat_cont::update_contacts update_nat_contacts;
  pauli::update_pairs update_pauli_pairs;
  pid::update_bundles update_pid_bundles;
  qa::update_free_pairs update_qa_pairs;

  pbar::render render_pbar;
  qa::prepare_nh prepare_nh;
  qa::process_candidates process_qa_candidates;

public:
  void thread_setup(thread_state &state);
};
} // namespace cg::simul