#include <cg/simul/thread.h>
namespace cg::simul {

void thread::main() {
  while (true) {
    step();
    if (st->cur_phase == phase::SIMUL_END)
      break;
  }
}

void thread::step() {
  switch (st->cur_phase) {
  case phase::SIMUL_INIT: {
#pragma omp single
    {
      st->simul_setup();
      st->traj_idx = 0;
      st->cur_phase = phase::TRAJ_INIT;
    }
    break;
  }
  case phase::TRAJ_INIT: {
    if (st->traj_idx < params->gen.num_of_traj) {
#pragma omp single
      {
        st->traj_equil_setup();
        st->until = min(params->gen.total_time, params->gen.equil_time);
        st->cur_phase = phase::EQUIL;
      }
      traj_equil_setup();
    } else {
#pragma omp single
      st->cur_phase = phase::SIMUL_END;
    }
    break;
  }
  case phase::EQUIL: {
    if (st->t < st->until) {
      advance_by_step();
    } else {
#pragma omp single
      {
        st->until = params->gen.total_time;
        st->cur_phase = phase::PROPER;
      }
    }
    break;
  }
  case phase::PROPER: {
    if (st->t < st->until) {
      advance_by_step();
    } else {
#pragma omp single
      {
        ++st->traj_idx;
        st->cur_phase = phase::TRAJ_INIT;
      }
    }
    break;
  }
  case phase::SIMUL_END:
    break;
  }
}

void thread::advance_by_step() {
  pre_eval_async();
  eval_forces();
  post_eval_async();
}

void thread::pre_eval_async() {
  st->dyn.omp_reset();

  if (params->qa.enabled) {
    prepare_nh.omp_async();
    if (st->ss_spec_crit)
      count_cys_neigh.omp_reset();
  }

  if (st->nl_required)
    nl_verify.omp_async();

#pragma omp barrier

  if (params->qa.enabled && st->ss_spec_crit)
    count_cys_neigh.omp_async();

  if (st->nl_invalid)
    fix_nl_async();

#pragma omp barrier
}

void thread::fix_nl_async() {
#pragma omp barrier

  switch (params->nl.algorithm) {
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
    if (st->pauli_enabled)
      update_pauli_pairs();
#pragma omp section
    if (st->nat_cont_enabled)
      update_nat_contacts();
#pragma omp section
    if (st->dh_enabled)
      update_dh_pairs();
#pragma omp section
    if (st->qa_enabled)
      update_qa_pairs();
#pragma omp section
    if (st->qa_enabled && st->ss_spec_crit)
      update_cys_neigh();
#pragma omp section
    if (st->pid_enabled)
      update_pid_bundles();
  }
}

void thread::eval_forces() {
  dyn.reset();

  if (st->chir_enabled)
    eval_chir_forces.omp_async();
  if (st->tether_enabled)
    eval_tether_forces.omp_async();
  if (st->lrep_enabled)
    eval_lrep_forces.omp_async();
  if (st->nat_cont_enabled)
    eval_nat_cont_forces.omp_async();
  if (st->pauli_enabled)
    eval_pauli_forces.omp_async();

  if (st->dh_enabled) {
    eval_const_dh_forces.omp_async();
    eval_rel_dh_forces.omp_async();
  }

  if (st->qa_enabled) {
    qa_loop_over_candidates.omp_async();
    process_qa_contacts.omp_async();
  }
  if (st->pid_enabled)
    eval_pid_forces.omp_async();

  if (st->nat_ang_enabled)
    eval_nat_ang_forces.omp_async();
  if (st->heur_ang_enabled)
    eval_heur_ang_forces.omp_async();
  if (st->heur_dih_enabled)
    eval_heur_dih_forces.omp_async();
  if (st->nat_dih_enabled) {
    eval_cnd_forces.omp_async();
    eval_snd_forces.omp_async();
  }

  if (st->solid_walls_enabled)
    eval_solid_wall_forces.omp_async();

  if (st->harmonic_walls_enabled) {
    hw_eval_free.omp_async();
    hw_eval_conn.omp_async();
  }

  if (st->lj_walls_enabled) {
    ljw_sift_free.omp_async();
    ljw_eval_conn.omp_async();
  }

  if (st->afm_enabled) {
    eval_vel_afm_forces.omp_async();
    eval_force_afm_forces.omp_async();
  }

  dyn.omp_reduce_v2(st->dyn, *this);
//  dyn.omp_reduce_v3(st->dyn, team);
#pragma omp barrier
}

void thread::post_eval_async() {
#pragma omp master
  {
    if (st->pbar_enabled)
      render_pbar();

    if (st->t > 0) {
      if (params->out.enabled)
        make_report();
    }

    if (st->dump_enabled)
      print_raw_data();

    if (st->ckpt_enabled)
      make_checkpoint();

    if (st->qa_enabled)
      qa_finish_processing();

    if (st->lj_walls_enabled)
      ljw_proc_cand();
  }

#pragma omp barrier

  if (st->lang_enabled) {
    switch (params->lang.type) {
    case lang::lang_type::NORMAL:
      lang_step.omp_async();
      break;
    case lang::lang_type::LEGACY:
      lang_legacy_step.omp_async();
      break;
    }
  }

  ++loop_idx;
}

} // namespace cg::simul