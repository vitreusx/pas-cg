#include <cg/simul/thread.h>
namespace cg::simul {

void thread::main() {
  do {
    loop();
  } while (st.is_running);
}

void thread::loop() {
  adjust_scenario();
  if (!st.is_running)
    return;

  pre_eval_async();
  eval_forces();
  post_eval_async();
}

void thread::adjust_scenario() {
#pragma omp barrier

#pragma omp master
  {
    if (!st.did_simul_setup) {
      st.simul_setup();
      st.did_simul_setup = true;
    }

    if (!st.did_traj_setup) {
      st.traj_setup();
      st.did_traj_setup = true;
    }

    if (st.t >= st.equil_time && !st.did_post_equil_setup) {
      st.post_equil_setup();
      st.did_post_equil_setup = true;
    }

    if (st.t >= st.total_time) {
      st.finish_trajectory();
      st.is_running = (st.traj_idx < params.gen.num_of_traj);
    }
  }
#pragma omp barrier

  if (!did_traj_setup) {
    traj_setup();
    did_traj_setup = true;
  }

  if (st.t >= st.total_time) {
    finish_trajectory();
    if (st.is_running) {
      adjust_scenario();
    } else {
      return;
    }
  }

  if (st.t >= st.equil_time && !did_post_equil_setup) {
    post_equil_setup();
    did_post_equil_setup = true;
  }
}

void thread::pre_eval_async() {
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
    fix_nl_async();

#pragma omp barrier
}

void thread::fix_nl_async() {
#pragma omp barrier

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
  if (params.lrep.enabled)
    eval_lrep_forces.omp_async();
  if (params.nat_cont.enabled)
    eval_nat_cont_forces.omp_async();

  if (params.nat_ang.enabled)
    eval_nat_ang_forces.omp_async();
  if (params.heur_ang.enabled)
    eval_heur_ang_forces.omp_async();

  if (params.heur_dih.enabled)
    eval_heur_dih_forces.omp_async();
  if (params.cnd.enabled)
    eval_cnd_forces.omp_async();
  if (params.snd.enabled)
    eval_snd_forces.omp_async();

  if (params.const_dh.enabled)
    eval_const_dh_forces.omp_async();
  if (params.rel_dh.enabled)
    eval_rel_dh_forces.omp_async();

  if (params.qa.enabled) {
    qa_loop_over_candidates.omp_async();
    process_qa_contacts.omp_async();
  }
  if (params.pid.enabled)
    eval_pid_forces.omp_async();

  if (params.pauli.enabled)
    eval_pauli_forces.omp_async();

  if (st.post_equil) {
    if (params.afm.enabled) {
      eval_vel_afm_forces.omp_async();
      eval_force_afm_forces.omp_async();
    }
  }

  dyn.omp_reduce_v2(st.dyn, *this);
#pragma omp barrier
}

void thread::post_eval_async() {
#pragma omp master
  {
    if (params.pbar.enabled)
      render_pbar();

    if (params.out.enabled)
      make_report();

    if (params.ckpt.enabled)
      make_checkpoint();

    if (params.qa.enabled)
      qa_finish_processing();
  }

#pragma omp barrier

  if (params.lang.enabled) {
    switch (params.lang.type) {
    case lang::lang_type::NORMAL:
      lang_step.omp_async();
      break;
    case lang::lang_type::LEGACY:
      lang_legacy_step.omp_async();
      break;
    }
  }

  //#ifdef COMPAT_MODE
  //#pragma omp master
  //  {
  //    if (params.out.enabled)
  //      make_report();
  //  };
  //#endif

  ++loop_idx;
}

} // namespace cg::simul