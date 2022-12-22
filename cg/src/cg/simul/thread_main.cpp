#include <cg/simul/thread.h>
namespace cg::simul {

void thread::main() {
  while (true) {
    step();
    if (st->cur_phase == phase::PROG_END)
      break;
  }
}

void thread::step() {
  if (st->trajectory_should_end()) {
#pragma omp barrier
#pragma omp single
    st->cur_phase = phase::TRAJ_END;
  }

  switch (st->cur_phase) {
  case phase::SIMUL_INIT:
    simul_init_step();
    break;
  case phase::TRAJ_INIT:
    traj_init_step();
    break;
  case phase::PULL_RELEASE:
    pull_release_step();
    break;
  case phase::EQUIL:
    equil_step();
    break;
  case phase::SQUEEZING:
    squeezing_step();
    break;
  case phase::REST_AFTER_SQUEEZING:
    rest_after_squeezing_step();
    break;
  case phase::FIND_FORCE_MIN:
    find_force_min_step();
    break;
  case phase::REST_AFTER_FORCE_MIN:
    rest_after_force_min_step();
    break;
  case phase::MAX_AMPLITUDE:
    max_amplitude_step();
    break;
  case phase::REST_AFTER_MAX_AMP:
    rest_after_max_amp_step();
    break;
  case phase::OSCILLATIONS:
    oscillations_step();
    break;
  case phase::REST_AFTER_OSCILLATIONS:
    rest_after_oscillations_step();
    break;
  case phase::FREEFORM:
    freeform_step();
    break;
  case phase::TRAJ_END:
    traj_end_step();
    break;
  case phase::SIMUL_END:
    simul_end_step();
    break;
  case phase::PROG_END:
    break;
  }
}

void thread::simul_init_step() {
#pragma omp barrier
#pragma omp single
  {
    st->simul_setup();
    st->traj_idx = 0;
    st->cur_phase = phase::TRAJ_INIT;

    setup_output();
    if (st->out_enabled)
      make_report.write_inputfile();
  }
}

void thread::traj_init_step() {
  if (st->traj_idx < params->gen.num_of_traj) {
#pragma omp barrier
#pragma omp single
    {
      st->traj_init();
      if (params->afm.type == "pull, then release") {
        st->until = min(params->gen.total_time, params->afm.pull_rel.time);
        st->cur_phase = phase::PULL_RELEASE;
        st->afm_enabled = true;
        st->setup_afm();
      } else {
        st->until = min(params->gen.total_time, params->gen.equil_time);
        st->cur_phase = phase::EQUIL;
      }
    }

    init_kernels();
  } else {
#pragma omp barrier
#pragma omp single
    st->cur_phase = phase::SIMUL_END;
  }
}

void thread::pull_release_step() {
  if (st->t < st->until) {
    advance_by_step();
  } else {
#pragma omp barrier
#pragma omp single
    {
      st->afm_enabled = false;
      st->until = st->t + params->gen.equil_time;
      st->cur_phase = phase::EQUIL;
    }
  }
}

void thread::equil_step() {
  if (st->t < st->until) {
    advance_by_step();
  } else {
#pragma omp barrier
#pragma omp single
    {
      if (params->afm.perform) {
        st->afm_enabled = true;
        st->setup_afm();
        st->cur_phase = phase::FREEFORM;
      } else {
        st->until = params->gen.total_time;
        st->since = st->t;
        st->cur_phase = phase::SQUEEZING;
      }
    }
    init_kernels();
  }
}

void thread::squeezing_step() {
#pragma omp barrier
#pragma omp single
  {
    if (!params->sbox.squeezing.perform) {
      st->since = st->t;
      st->cur_phase = phase::FIND_FORCE_MIN;
    } else {
      auto const &sp = params->sbox.squeezing;
      auto ext = st->box.extent();
      auto cur_vol = ext.x() * ext.y() * ext.z();
      auto target_vol = st->num_res / (real)sp.target_density;

      if (cur_vol <= target_vol) {
        st->until = min(st->t + (real)params->sbox.rest_period,
                        (real)params->gen.total_time);
        st->cur_phase = phase::REST_AFTER_SQUEEZING;
      } else {
        real vel;
        if (cur_vol > (real)2.0 * target_vol) {
          real accel_time = params->sbox.accel_dist / sp.vel_above_2V;
          real time_frac = min((st->t - st->since) / accel_time, (real)1.0);
          vel = time_frac * sp.vel_above_2V;
        } else {
          real accel_time = params->sbox.accel_dist / params->sbox.target_vel;
          real time_frac = min((st->t - st->since) / accel_time, (real)1.0);
          vel = time_frac * params->sbox.target_vel;
        }

        real shift = (real)0.5 * vel * params->lang.dt;
        st->adjust_wall_pos(vec3r(-shift, -shift, -shift), vec3r::Zero());
      }
    }
  }

  if (st->cur_phase == phase::SQUEEZING)
    advance_by_step();
}

void thread::rest_after_squeezing_step() {
  if (st->t < st->until) {
    advance_by_step();
  } else {
#pragma omp barrier
#pragma omp single
    st->cur_phase = phase::FIND_FORCE_MIN;
  }
}

void thread::find_force_min_step() {
#pragma omp barrier
#pragma omp single
  {
    if (!params->sbox.force_min.perform) {
      st->cur_phase = phase::FREEFORM;
    } else {
      if (st->avg_z_force->has_value() &&
          abs(st->avg_z_force->value()) < 5e-2 + 5e-5 * st->num_res) {
        st->until = st->t + params->sbox.rest_period;
        st->cur_phase = phase::REST_AFTER_FORCE_MIN;
      } else {
        auto const &fmp = params->sbox.force_min;
        real accel_time = params->sbox.accel_dist / params->sbox.target_vel;
        real time_vel_frac = min((st->t - st->since) / accel_time, (real)1.0);
        real force_vel_frac =
            st->avg_z_force->value_or(0) / fmp.force_for_max_vel;
        force_vel_frac =
            copysign(min(abs(force_vel_frac), (real)1.0), force_vel_frac);
        real vel = params->sbox.target_vel * time_vel_frac * force_vel_frac;

        real shift = (real)0.5 * vel * params->lang.dt;
        st->adjust_wall_pos(vec3r(0, 0, shift), vec3r::Zero());
      }
    }
  }

  if (st->cur_phase == phase::FIND_FORCE_MIN)
    advance_by_step();
}

void thread::rest_after_force_min_step() {
  if (st->t < st->until) {
    advance_by_step();
  } else {
#pragma omp barrier
#pragma omp single
    {
      auto const &op = params->sbox.oscillations;
      if (op.amplitude.variant == "absolute") {
        st->amplitude = op.amplitude.abs_value;
      } else {
        st->shear = op.type == "shear";
        real base_val = st->shear ? st->box.extent().x() : st->box.extent().z();
        st->amplitude = op.amplitude.rel_value * base_val;
      }
      st->since = st->t;
      st->cur_phase = phase::MAX_AMPLITUDE;
    }
  }
}

void thread::max_amplitude_step() {
#pragma omp barrier
#pragma omp single
  {
    if (!params->sbox.oscillations.perform) {
      st->cur_phase = phase::FREEFORM;
    } else {
      real accel_time = params->sbox.accel_dist / params->sbox.target_vel;
      real time_frac = min((st->t - st->since) / accel_time, (real)1.0);
      real vel = time_frac * params->sbox.target_vel;
      real shift = vel * params->lang.dt;

      if (st->displacement > (real)0.5 * st->amplitude) {
        st->until = st->t + params->sbox.rest_period;
        st->cur_phase = phase::REST_AFTER_MAX_AMP;
      } else {
        st->displacement += shift;
        if (st->shear) {
          if (st->walls[NEG_Z])
            st->walls[NEG_Z]->plane.normal().x() += shift;
          if (st->walls[POS_Z])
            st->walls[POS_Z]->plane.normal().x() += shift;
        } else {
          st->adjust_wall_pos(vec3r(0, 0, -shift), vec3r::Zero());
        }
      }
    }
  }

  if (st->cur_phase == phase::MAX_AMPLITUDE)
    advance_by_step();
}

void thread::rest_after_max_amp_step() {
  if (st->t < st->until) {
    advance_by_step();
  } else {
#pragma omp barrier
#pragma omp single
    {
      st->since = st->t;
      auto const &op = params->sbox.oscillations;
      st->omega = op.angular_freq;
      st->until =
          st->t + ((real)M_PI * 2 * op.num_of_cycles + M_PI_2) / st->omega;
      st->cur_phase = phase::OSCILLATIONS;
    }
  }
}

void thread::oscillations_step() {
#pragma omp barrier
#pragma omp single
  {
    if (st->t >= st->until) {
      st->until = st->t + params->sbox.rest_period;
      st->cur_phase = phase::REST_AFTER_OSCILLATIONS;
    } else {
      real vel = (real)0.5 * st->amplitude * st->omega *
                 sin(st->omega * (st->t - st->since));
      real shift = vel * params->lang.dt;

      if (st->shear) {
        if (st->walls[NEG_Z])
          st->walls[NEG_Z]->plane.normal().x() += shift;
        if (st->walls[POS_Z])
          st->walls[POS_Z]->plane.normal().x() += shift;
      } else {
        st->adjust_wall_pos(vec3r(0, 0, shift), vec3r::Zero());
      }
    }
  }

  if (st->cur_phase == phase::OSCILLATIONS)
    advance_by_step();
}

void thread::rest_after_oscillations_step() {
  if (st->t < st->until) {
    advance_by_step();
  } else {
#pragma omp barrier
#pragma omp single
    st->cur_phase = phase::FREEFORM;
  }
}

void thread::freeform_step() {
  advance_by_step();
}

void thread::traj_end_step() {
#pragma omp barrier
#pragma omp single
  {
    if (params->nat_cont.unfolding_study.measure_times) {
      auto num_nc = st->all_native_contacts.size();
      auto &traj_nc_times = st->nc_times.emplace_back(num_nc);
      for (int idx = 0; idx < num_nc; ++idx)
        traj_nc_times[idx] = st->all_native_contacts[idx].change_t();

      if (st->num_changed == num_nc) {
        real unfold_t = 0.0;
        for (auto const &cont : st->all_native_contacts)
          unfold_t = max(unfold_t, cont.change_t());
        st->nc_unfold_times.push_back(unfold_t);
      }
    }

    if (st->out_enabled)
      make_report.emit_all();

    ++st->traj_idx;
    st->cur_phase = phase::TRAJ_INIT;
  }
}

void thread::simul_end_step() {
#pragma omp barrier
#pragma omp single
  {
    if (st->out_enabled)
      make_report.at_simul_end();

    st->cur_phase = phase::PROG_END;
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
#pragma omp single nowait
      reset_cys_neigh();
  }

  if (st->def_nl.required)
    def_nl.verify.omp_async();

  if (st->long_range_nl.required)
    long_range_nl.verify.omp_async();

#pragma omp barrier

  if (params->qa.enabled && st->ss_spec_crit)
    count_cys_neigh.omp_async();

  if (st->def_nl.invalid || st->long_range_nl.invalid) {
    if (st->def_nl.invalid)
      fix_def_nl_async();
    if (st->long_range_nl.invalid)
      fix_lr_nl_async();
#pragma omp barrier
  }
}

void thread::fix_def_nl_async() {
  // #pragma omp barrier
  //   def_nl.legacy.omp_async();

  def_nl.cell.omp_async();
#pragma omp barrier

#pragma omp master
  {
    if (st->nat_cont_enabled)
      update_nat_contacts();
    if (st->qa_enabled && st->ss_spec_crit)
      update_cys_neigh();
    if (st->qa_enabled)
      update_qa_pairs();
    if (st->pid_enabled)
      update_pid_bundles();
    if (st->pauli_enabled)
      update_pauli_pairs();
  };

#pragma omp barrier
  eval_divs.update();
}

void thread::fix_lr_nl_async() {
  // #pragma omp barrier
  //   long_range_nl.legacy.omp_async();

  long_range_nl.cell.omp_async();
#pragma omp barrier

#pragma omp master
  {
    auto nl_ = &st->long_range_nl.nl;
    if (st->nat_cont_enabled)
      update_nat_contacts.mark_as_taken(nl_);
    if (st->qa_enabled)
      update_qa_pairs.mark_as_taken(nl_);
    if (st->dh_enabled)
      update_dh_pairs();
  };

#pragma omp barrier
  eval_divs.update();
}

template <typename E>
std::ostream &operator<<(std::ostream &os, vec3_expr<E> const &e) {
  os << "(" << e.x() << ", " << e.y() << ", " << e.z() << ")";
  return os;
}

void thread::eval_forces() {
  dyn.reset();

  //  if (st->chir_enabled)
  //    eval_chir_forces.omp_async();
  //  if (st->tether_enabled)
  //    eval_tether_forces.omp_async();
  //  if (st->lrep_enabled)
  //    eval_lrep_forces.omp_async();
  //  if (st->nat_cont_enabled)
  //    eval_nat_cont_forces.omp_async();
  //  if (st->pauli_enabled)
  //    eval_pauli_forces.omp_async();
  //
  //  if (st->dh_enabled) {
  //    eval_const_dh_forces.omp_async();
  //    eval_rel_dh_forces.omp_async();
  //  }
  //
  //  if (st->qa_enabled) {
  //    qa_loop_over_candidates.omp_async();
  //    process_qa_contacts.omp_async();
  //  }
  //  if (st->pid_enabled)
  //    eval_pid_forces.omp_async();
  //
  //  if (st->nat_ang_enabled)
  //    eval_nat_ang_forces.omp_async();
  //  if (st->heur_ang_enabled)
  //    eval_heur_ang_forces.omp_async();
  //  if (st->heur_dih_enabled)
  //    eval_heur_dih_forces.omp_async();
  //  if (st->nat_dih_enabled) {
  //    eval_cnd_forces.omp_async();
  //    eval_snd_forces.omp_async();
  //  }
  //
  //  if (st->solid_walls_enabled)
  //    eval_solid_wall_forces.omp_async();
  //
  //  if (st->harmonic_walls_enabled) {
  //    hw_eval_free.omp_async();
  //    hw_eval_conn.omp_async();
  //  }
  //
  //  if (st->lj_walls_enabled) {
  //    ljw_sift_free.omp_async();
  //    ljw_eval_conn.omp_async();
  //  }
  //
  //  if (st->afm_enabled) {
  //    eval_vel_afm_forces.omp_async();
  //    eval_force_afm_forces.omp_async();
  //  }

  eval_divs.omp_async();

  dyn.omp_reduce_v2(st->dyn, *this);
  //  dyn.omp_reduce_v3(st->dyn, *team);
//  dyn.omp_reduce_v4(st->dyn, team->v4_shared_, dyn_v4_priv_);
#pragma omp barrier
}

void thread::post_eval_async() {
#pragma omp master
  {
    if (st->pbar_enabled)
      render_pbar();

    if (st->out_enabled && st->t > 0)
      make_report();

    if (st->dump_enabled)
      print_raw_data();

    if (st->ckpt_enabled)
      make_checkpoint();

    if (st->qa_enabled)
      qa_finish_processing();

    if (st->lj_walls_enabled)
      ljw_proc_cand();

    if (log_wall_forces_enabled)
      log_wall_forces();
  }

#pragma omp barrier

  if (st->lang_enabled) {
    //    lang_step.omp_async();
    lang_fast_step.omp_async();
  }

  ++loop_idx;
}

} // namespace cg::simul