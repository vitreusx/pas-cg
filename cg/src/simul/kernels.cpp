#include "simul/kernels.h"
using namespace cg::simul;

void kernels::thread_setup(thread_state &state) {
  eval_chir_forces.F = state.dyn.F.view();
  eval_chir_forces.V = &state.dyn.V;
  eval_const_dh_forces.F = state.dyn.F.view();
  eval_const_dh_forces.V = &state.dyn.V;
  eval_rel_dh_forces.F = state.dyn.F.view();
  eval_rel_dh_forces.V = &state.dyn.V;
  eval_fafm_forces.F = state.dyn.F.view();
  eval_heur_ang_forces.F = state.dyn.F.view();
  eval_heur_ang_forces.V = &state.dyn.V;
  eval_heur_dih_forces.F = state.dyn.F.view();
  eval_heur_dih_forces.V = &state.dyn.V;
  eval_nat_ang_forces.F = state.dyn.F.view();
  eval_nat_ang_forces.V = &state.dyn.V;
  eval_nat_cont_forces.F = state.dyn.F.view();
  eval_nat_cont_forces.V = &state.dyn.V;
  eval_cnd_forces.F = state.dyn.F.view();
  eval_cnd_forces.V = &state.dyn.V;
  eval_snd_forces.F = state.dyn.F.view();
  eval_snd_forces.V = &state.dyn.V;
  eval_pauli_forces.F = state.dyn.F.view();
  eval_pauli_forces.V = &state.dyn.V;
  eval_pid_forces.F = state.dyn.F.view();
  eval_pid_forces.V = &state.dyn.V;
  process_qa_contacts.F = state.dyn.F.view();
  process_qa_contacts.V = &state.dyn.V;
  eval_tether_forces.F = state.dyn.F.view();
  eval_tether_forces.V = &state.dyn.V;
  eval_vafm_forces.F = state.dyn.F.view();
  lang_step.gen = &state.gen;
}