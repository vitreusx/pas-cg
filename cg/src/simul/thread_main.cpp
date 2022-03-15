#include "simul/thread.h"
namespace cg::simul {

void thread::main() {
  do {
    adjust_scenario();
    if (!st.is_running)
      break;

    pre_eval();
    eval_forces();
    post_eval();
  } while (true);
}

void thread::adjust_scenario() {
  if (!did_overall_setup) {
    if (!st.did_overall_setup) {
#pragma omp master
      {
        st.overall_setup();
        st.did_overall_setup = true;
      }
#pragma omp barrier
    }

    overall_setup();
    did_overall_setup = true;
  }

  if (st.t >= st.total_time) {
    st.is_running = false;
    return;
  }

  if (st.t >= st.equil_time && !did_post_equil_setup) {
    if (!st.did_post_equil_setup) {
#pragma omp master
      {
        st.post_equil_setup();
        st.did_post_equil_setup = true;
      }
#pragma omp barrier
    }

    post_equil_setup();
    did_post_equil_setup = true;
  }

#pragma omp barrier
}

void thread::pre_eval() {
  st.dyn.omp_reset();

  if (params.qa.enabled) {
    prepare_nh.omp_async();
    if (st.ss_spec_crit)
      count_cys_neigh.omp_reset();
  }

  if (st.nl_required)
    nl_verify.omp_async();

#pragma omp barrier
  if (params.qa.enabled && st.ss_spec_crit)
    count_cys_neigh.omp_async();

  if (st.nl_invalid)
    fix_nl();

#pragma omp barrier
}

void thread::fix_nl() {
  switch (params.nl.algorithm) {
  case nl::parameters::CELL:
#pragma omp master
    nl_cell();
    break;
  case nl::parameters::LEGACY:
    nl_legacy.omp_async();
    break;
  }
#pragma omp barrier

#pragma omp sections
  {
#pragma omp section
    if (params.pauli.enabled)
      update_pauli_pairs();
#pragma omp section
    if (params.nat_cont.enabled)
      update_nat_contacts();
#pragma omp section
    if (params.const_dh.enabled || params.rel_dh.enabled)
      update_dh_pairs();
#pragma omp section
    if (params.qa.enabled)
      update_qa_pairs();
#pragma omp section
    if (params.qa.enabled && st.ss_spec_crit)
      update_cys_neigh();
#pragma omp section
    if (params.pid.enabled)
      update_pid_bundles();
  }
}

void thread::eval_forces() {
  dyn.reset();

  if (params.chir.enabled)
    eval_chir_forces.omp_async();
  if (params.tether.enabled)
    eval_tether_forces.omp_async();
  if (params.nat_ang.enabled)
    eval_nat_ang_forces.omp_async();
  if (params.cnd.enabled)
    eval_cnd_forces.omp_async();
  if (params.snd.enabled)
    eval_snd_forces.omp_async();
  if (params.pauli.enabled)
    eval_pauli_forces.omp_async();
  if (params.nat_cont.enabled)
    eval_nat_cont_forces.omp_async();
  if (params.const_dh.enabled)
    eval_const_dh_forces.omp_async();
  if (params.rel_dh.enabled)
    eval_rel_dh_forces.omp_async();
  if (params.heur_ang.enabled)
    eval_heur_ang_forces.omp_async();
  if (params.heur_dih.enabled)
    eval_heur_dih_forces.omp_async();
  if (params.qa.enabled) {
    sift_qa_candidates.omp_async();
    process_qa_contacts.omp_async();
  }
  if (params.pid.enabled)
    eval_pid_forces.omp_async();

  if (st.post_equil) {
    if (params.afm.enabled) {
      eval_vel_afm_forces.omp_async();
      eval_force_afm_forces.omp_async();
    }
  }

  dyn.omp_reduce(st.dyn);
#pragma omp barrier
}

void thread::post_eval() {
#pragma omp master
  {
    if (params.pbar.enabled) {
      render_pbar();
    }

    if (params.out.enabled) {
      make_report();
    }

    if (params.qa.enabled) {
      qa_finish_processing();
    }
  }

  if (params.lang.enabled) {
#pragma omp barrier
    switch (params.lang.type) {
    case lang::lang_type::NORMAL:
      lang_step.omp_async();
      break;
    case lang::lang_type::LEGACY:
      lang_legacy_step.omp_async();
      break;
    }
#pragma omp barrier
  }
}
} // namespace cg::simul