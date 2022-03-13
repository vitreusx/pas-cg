#include "simul/thread.h"
using namespace cg::simul;

void thread::overall_setup() {
  setup_gen();
  setup_dyn();
  setup_output();
  setup_langevin();
  setup_pbar();
  setup_nl();
  setup_chir();
  setup_tether();
  setup_nat_ang();
  setup_nat_dih();
  setup_heur_dih();
  setup_pauli();
  setup_nat_cont();
  setup_dh();
  setup_qa();
  setup_pid();
}

void thread::setup_gen() {
  if (params.gen.debug_mode) {
    feclearexcept(FE_ALL_EXCEPT);
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
  }
}

void thread::setup_dyn() { dyn = dynamics(st.num_res); }

void thread::setup_output() {
  if (params.out.enabled) {
    make_report.period = params.out.period;
    make_report.t = &st.t;
    make_report.last_t = &st.report_last_t;
    make_report.hooks = &hooks;
    make_report.state = &st.report;

    export_pdb.ref_model = &st.model;
    export_pdb.r = st.r.view();
    export_pdb.res_map = &st.res_map;
    hooks.emplace_back(&export_pdb);

    add_stats.t = &st.t;
    add_stats.V = &st.dyn.V;
    add_stats.atype = st.atype.view();
    add_stats.mass = st.comp_aa_data.mass.view();
    add_stats.chain_first = st.chain_first.view();
    add_stats.chain_last = st.chain_last.view();
    add_stats.r = st.r.view();
    hooks.emplace_back(&add_stats);

    add_structure.res_map = &st.res_map;
    add_structure.r = st.r.view();
    add_structure.model = &st.model;
    hooks.emplace_back(&add_structure);

    report_gyr.r = st.r.view();
    hooks.emplace_back(&report_gyr);

    compute_rmsd.orig_r = st.orig_r.view();
  }
}

void thread::setup_langevin() {

  if (params.lang.enabled) {
    auto &step = params.lang.type == lang::lang_type::NORMAL
                     ? (lang::step_base &)lang_step
                     : (lang::step_base &)lang_legacy_step;

    step.temperature = params.lang.temperature;
    step.t = &st.t;
    step.dt = params.lang.dt;
    step.gamma_factor = params.lang.gamma;
    step.gen = &st.gen;
    step.mass = st.comp_aa_data.mass.view();
    step.num_particles = st.num_res;
    step.atype = st.atype.view();
    step.r = st.r.view();
    step.mass_inv = st.mass_inv.view();
    step.mass_rsqrt = st.mass_rsqrt.view();

    step.v = st.v.view();
    if (params.out.enabled)
      add_stats.v = st.v.view();

    step.y0 = st.y0.view();
    step.y1 = st.y1.view();
    step.y2 = st.y2.view();
    step.y3 = st.y3.view();
    step.y4 = st.y4.view();
    step.y5 = st.y5.view();

    step.true_t = &st.true_t;

    step.F = st.dyn.F.view();
  }
}

void thread::setup_pbar() {
  if (params.pbar.enabled) {
    auto &render = render_pbar;
    render.width = params.pbar.width;
    render.total_time = params.gen.total_time;
    render.period_s = params.pbar.update_period.in("s");
    render.t = &st.t;
    render.V = &st.dyn.V;
    render.start_clock = &st.pbar_start_clock;
    render.last_clock = &st.pbar_last_clock;
    render.t = &st.t;

    render.is_first = &st.pbar_first_time;
  }
}

void thread::setup_nl() {
  if (params.nl.algorithm == nl::parameters::LEGACY) {
    auto &legacy = nl_legacy;
    legacy.pad_factor = params.nl.pad_factor;
    legacy.r = st.r.view();
    legacy.simul_box = &st.box;
    legacy.t = &st.t;
    legacy.chain_idx = st.chain_idx.view();
    legacy.seq_idx = st.seq_idx.view();
    legacy.num_particles = st.num_res;
    legacy.nl_data = &st.nl;
    legacy.invalid = &st.nl_invalid;
    legacy.max_cutoff = &st.max_cutoff;
  } else {
    auto &cell = nl_cell;
    cell.pad_factor = params.nl.pad_factor;
    cell.r = st.r.view();
    cell.simul_box = &st.box;
    cell.t = &st.t;
    cell.chain_idx = st.chain_idx.view();
    cell.seq_idx = st.seq_idx.view();
    cell.num_particles = st.num_res;
    cell.nl_data = &st.nl;
    cell.invalid = &st.nl_invalid;
    cell.max_cutoff = &st.max_cutoff;
    cell.res_cell_idx = st.res_cell_idx.view();
    cell.reordered_idx = st.reordered_idx.view();
    cell.num_res_in_cell = &st.num_res_in_cell;
    cell.cell_offset = &st.cell_offset;
    cell.all_pairs = &st.all_pairs;
  }

  auto &verify = nl_verify;
  verify.r = st.r.view();
  verify.nl_data = &st.nl;
  verify.simul_box = &st.box;
  verify.invalid = &st.nl_invalid;
  verify.first_time = &st.verify_first_time;
  verify.num_particles = st.num_res;
}

void thread::setup_chir() {
  if (params.chir.enabled) {
    auto &eval = eval_chir_forces;
    eval.e_chi = params.chir.e_chi;
    eval.quads = st.chir_quads.view();
    eval.V = &dyn.V;
    eval.F = dyn.F.view();
  }
}

void thread::setup_tether() {
  if (params.tether.enabled) {
    auto &eval = eval_tether_forces;
    eval.H1 = params.tether.H1;
    eval.H2 = params.tether.H2;
    eval.def_length = params.tether.def_length;
    eval.r = st.r.view();
    eval.tethers = st.tether_pairs.view();
    eval.V = &dyn.V;
    eval.F = dyn.F.view();
  }
}

void thread::setup_nat_ang() {
  if (params.nat_ang.enabled) {
    auto &eval = eval_nat_ang_forces;
    eval.k = params.nat_ang.k;
    eval.r = st.r.view();
    eval.angles = st.native_angles.view();
    eval.V = &dyn.V;
    eval.F = dyn.F.view();
  }
}

void thread::setup_heur_ang() {
  if (params.heur_ang.enabled) {
    auto &eval = eval_heur_ang_forces;
    eval.r = st.r.view();
    eval.angles = st.heur_angles.view();
    eval.V = &dyn.V;
    eval.F = dyn.F.view();

    for (auto const &heur_pair : aa_heur_pair::all()) {
      for (int d = 0; d <= heur_ang::eval_forces::POLY_DEG; ++d) {
        eval.poly_coeffs[d][(uint8_t)heur_pair] =
            params.heur_ang.coeffs.at(heur_pair).poly[d];
      }
    }
  }
}

void thread::setup_nat_dih() {
  if (params.cnd.enabled) {
    auto &eval = eval_cnd_forces;
    eval.CDA = params.cnd.CDA;
    eval.CDB = params.cnd.CDB;
    eval.r = st.r.view();
    eval.dihedrals = st.native_dihedrals.view();
    eval.V = &dyn.V;
    eval.F = dyn.F.view();
  }

  if (params.snd.enabled) {
    auto &eval = eval_snd_forces;
    eval.CDH = params.snd.CDH;
    eval.r = st.r.view();
    eval.dihedrals = st.native_dihedrals.view();
    eval.V = &dyn.V;
    eval.F = dyn.F.view();
  }
}

void thread::setup_heur_dih() {
  if (params.heur_dih.enabled) {
    auto &eval = eval_heur_dih_forces;
    eval.r = st.r.view();
    eval.V = &dyn.V;
    eval.F = dyn.F.view();

    for (auto const &heur_pair : aa_heur_pair::all()) {
      auto idx = (uint8_t)heur_pair;
      auto const &coeffs = params.heur_dih.coeffs.at(heur_pair);
      eval.coeffs.const_[idx] = coeffs.const_;
      eval.coeffs.sin[idx] = coeffs.sin;
      eval.coeffs.cos[idx] = coeffs.cos;
      eval.coeffs.sin2[idx] = coeffs.sin2;
      eval.coeffs.cos2[idx] = coeffs.cos2;
      eval.coeffs.sin_cos[idx] = coeffs.sin_cos;
    }
  }
}

void thread::setup_pauli() {
  if (params.pauli.enabled) {
    auto &eval = eval_pauli_forces;
    eval.r_excl = params.pauli.r_excl;
    eval.depth = params.pauli.depth;
    eval.r = st.r.view();
    eval.simul_box = &st.box;
    eval.pairs = &st.pauli_pairs;
    eval.V = &dyn.V;
    eval.F = dyn.F.view();

    auto &update = update_pauli_pairs;
    update.r = st.r.view();
    update.simul_box = &st.box;
    update.nl = &st.nl;
    update.pairs = &st.pauli_pairs;
    update.r_excl = st.pauli_cutoff;
  }
}

void thread::setup_nat_cont() {
  if (params.nat_cont.enabled) {
    nl_legacy.all_nat_cont = st.nat_cont_excl.view();
    nl_cell.all_nat_cont = st.nat_cont_excl.view();

    auto &eval = eval_nat_cont_forces;
    eval.depth = params.nat_cont.lj_depth;
    eval.simul_box = &st.box;
    eval.contacts = &st.cur_native_contacts;
    eval.r = st.r.view();
    eval.active_thr = params.nat_cont.active_thr;
    eval.t = &st.t;
    eval.all_contacts = st.all_native_contacts.view();
    eval.V = &dyn.V;
    eval.F = dyn.F.view();

    auto &update_ref = update_nat_contacts;
    update_ref.r = st.r.view();
    update_ref.simul_box = &st.box;
    update_ref.nl = &st.nl;
    update_ref.all_contacts = st.all_native_contacts.view();
    update_ref.contacts = &st.cur_native_contacts;

    if (params.out.enabled) {
      auto &report_nc = report_nc_stuff;
      report_nc.all_contacts = st.all_native_contacts.view();
      report_nc.params = &params.nat_cont;
      report_nc.r = st.r.view();
      report_nc.chain_idx = st.chain_idx.view();
      hooks.push_back(&report_nc);
    }
  }
}

void thread::setup_dh() {
  auto &update = update_dh_pairs;
  update.r = st.r.view();
  update.simul_box = &st.box;
  update.nl = &st.nl;
  update.pairs = &st.dh_pairs;
  update.atype = st.atype.view();
  for (auto const &aa : amino_acid::all())
    update.q[(uint8_t)aa] = st.comp_aa_data.charge[(uint8_t)aa];

  if (params.const_dh.enabled) {
    auto &eval = eval_const_dh_forces;
    eval.set_V_factor(params.const_dh.permittivity);
    eval.screen_dist_inv = 1.0 / params.const_dh.screening_dist;
    eval.r = st.r.view();
    eval.simul_box = &st.box;
    eval.es_pairs = &st.dh_pairs;
    eval.V = &dyn.V;
    eval.F = dyn.F.view();
    update.cutoff = 2.0 * params.const_dh.screening_dist;
  }

  if (params.rel_dh.enabled) {
    auto &eval = eval_rel_dh_forces;
    eval.set_V_factor(params.rel_dh.perm_factor);
    eval.screen_dist_inv = 1.0 / params.rel_dh.screening_dist;
    eval.r = st.r.view();
    eval.simul_box = &st.box;
    eval.es_pairs = &st.dh_pairs;
    eval.V = &dyn.V;
    eval.F = dyn.F.view();
    update.cutoff = 2.0 * params.rel_dh.screening_dist;
  }
}

void thread::setup_qa() {
  if (params.qa.enabled) {
    auto &prep_nh = prepare_nh;
    prep_nh.r = st.r.view();
    prep_nh.num_particles = st.num_res;
    prep_nh.n = st.n.view();
    prep_nh.h = st.h.view();
    prep_nh.prev = st.prev.view();
    prep_nh.next = st.next.view();

    auto &sift = sift_qa_candidates;
    sift.r = st.r.view();
    sift.atype = st.atype.view();
    sift.simul_box = &st.box;
    sift.free_pairs = &st.qa_free_pairs;
    sift.sync = st.sync_values.view();
    sift.candidates = &st.qa_candidates;
    sift.n = st.n.view();
    sift.h = st.h.view();

    sift.min_abs_cos_hr = params.qa.min_cos_hr;
    sift.min_abs_cos_hh = params.qa.min_cos_hh;
    sift.max_cos_nr = params.qa.max_cos_nr;
    sift.formation_tolerance = params.qa.formation_tolerance;

    sift.req_min_dist[(int16_t)qa::contact_type::BACK_BACK()] =
        params.lj_vars.bb.r_min();
    sift.req_min_dist[(int16_t)qa::contact_type::BACK_SIDE()] =
        params.lj_vars.bs.r_min();
    sift.req_min_dist[(int16_t)qa::contact_type::SIDE_BACK()] =
        params.lj_vars.sb.r_min();

    for (auto const &aa1 : amino_acid::all()) {
      for (auto const &aa2 : amino_acid::all()) {
        sift.req_min_dist[(int16_t)qa::contact_type::SIDE_SIDE(aa1, aa2)] =
            params.lj_vars.ss.at({aa1, aa2}).r_max();
      }
    }

    for (auto const &aa : amino_acid::all())
      sift.ptype[(uint8_t)aa] = st.comp_aa_data.ptype[(uint8_t)aa];

    auto &fin_proc = qa_finish_processing;
    fin_proc.candidates = &st.qa_candidates;
    fin_proc.contacts = &st.qa_contacts;
    fin_proc.sync = st.sync_values.view();
    fin_proc.free_pairs = &st.qa_free_pairs;
    fin_proc.t = &st.t;
    fin_proc.removed = &st.qa_removed;
    fin_proc.num_contacts = &st.num_qa_contacts;

    auto &proc_cont = process_qa_contacts;
    proc_cont.cycle_time = params.qa.phase_dur;
    proc_cont.cycle_time_inv = 1.0 / proc_cont.cycle_time;
    proc_cont.ljs = qa::lj_variants(params.lj_vars);
    proc_cont.set_factor(params.qa.breaking_factor);
    proc_cont.t = &st.t;
    proc_cont.sync = st.sync_values.view();
    proc_cont.contacts = &st.qa_contacts;
    proc_cont.simul_box = &st.box;
    proc_cont.r = st.r.view();
    proc_cont.free_pairs = &st.qa_free_pairs;
    proc_cont.removed = &st.qa_removed;
    proc_cont.V = &dyn.V;
    proc_cont.F = dyn.F.view();

    auto &update = update_qa_pairs;
    update.r = st.r.view();
    update.simul_box = &st.box;
    update.nl = &st.nl;
    update.pairs = &st.qa_free_pairs;

    if (params.qa.disulfide.has_value()) {
      if (st.ss_spec_crit) {
        auto const &spec_crit_params = params.qa.disulfide->spec_crit;

        sift.part_of_ssbond = st.part_of_ssbond.view();
        sift.disulfide_special_criteria = st.ss_spec_crit;
        sift.neigh = st.neigh_count.view();

        fin_proc.disulfide_special_criteria = st.ss_spec_crit;
        fin_proc.part_of_ssbond = st.part_of_ssbond.view();

        proc_cont.disulfide_special_criteria = st.ss_spec_crit;
        proc_cont.disulfide = params.qa.disulfide->force;

        proc_cont.neigh = st.neigh_count.view();
        proc_cont.max_neigh_count = spec_crit_params.max_neigh_count;
        proc_cont.part_of_ssbond = st.part_of_ssbond.view();

        auto &update_cys = update_cys_neigh;
        update_cys.neigh = &st.qa_cys_neigh;
        update_cys.neigh_radius = spec_crit_params.neigh_radius;
        update_cys.nl = &st.nl;
        update_cys.simul_box = &st.box;
        update_cys.r = st.r.view();

        auto &count_cys = count_cys_neigh;
        count_cys.simul_box = &st.box;
        count_cys.neigh_radius = spec_crit_params.neigh_radius;
        count_cys.neigh = &st.qa_cys_neigh;
        count_cys.neigh_count = st.neigh_count.view();
        count_cys.cys_indices = st.cys_indices.view();
        count_cys.r = st.r.view();

#pragma omp single nowait
        st.max_cutoff = max(st.max_cutoff, (real)spec_crit_params.neigh_radius);
      }
    }

    real cutoff = 0.0;
    for (int idx = 0; idx < qa::contact_type::NUM_TYPES; ++idx)
      cutoff = max(cutoff, sift.req_min_dist[idx]);

#pragma omp single nowait
    st.max_cutoff = max(st.max_cutoff, cutoff);

    if (params.out.enabled) {
      auto &report_qa = report_qa_stuff;
      report_qa.dyn_ss = params.qa.disulfide.has_value();
      report_qa.sync_values = st.sync_values.view();
      report_qa.contacts = &st.qa_contacts;
      report_qa.process_cont = &process_qa_contacts;
      report_qa.chain_idx = st.chain_idx.view();
      hooks.push_back(&report_qa);
    }
  }
}

void thread::setup_pid() {
  if (params.pid.enabled) {
    auto &eval = eval_pid_forces;
    eval.V = &dyn.V;
    eval.F = dyn.F.view();

    eval.bb_plus_lam.version() = params.pid.lambda_version;
    eval.bb_plus_lam.alpha() = params.pid.bb_plus.alpha;
    eval.bb_plus_lam.psi_0() = params.pid.bb_plus.psi_0;
    eval.bb_plus_lj.r_min() = params.pid.bb_plus.r_min;
    eval.bb_plus_lj.depth() = params.pid.bb_plus.depth;

    eval.bb_minus_lam.version() = params.pid.lambda_version;
    eval.bb_minus_lam.alpha() = params.pid.bb_minus.alpha;
    eval.bb_minus_lam.psi_0() = params.pid.bb_minus.psi_0;
    eval.bb_minus_lj.r_min() = params.pid.bb_minus.r_min;
    eval.bb_minus_lj.depth() = params.pid.bb_minus.depth;

    eval.ss_lam.alpha() = params.pid.ss.alpha;
    eval.ss_lam.psi_0() = params.pid.ss.psi_0;
    eval.ss_ljs = st.ss_ljs.view();

    eval.r = st.r.view();
    eval.simul_box = &st.box;
    eval.bundles = &st.pid_bundles;
    eval.prev = st.prev.view();
    eval.next = st.next.view();
    eval.V = &dyn.V;
    eval.F = dyn.F.view();

    auto &update = update_pid_bundles;
    update.r = st.r.view();
    update.prev = st.prev.view();
    update.next = st.next.view();
    update.atype = st.atype.view();
    update.simul_box = &st.box;
    update.nl = &st.nl;
    update.bundles = &st.pid_bundles;

    real cutoff = 0.0;
    cutoff = max(cutoff, params.lj_vars.bb.cutoff());
    cutoff = max(cutoff, params.lj_vars.bs.cutoff());
    cutoff = max(cutoff, params.lj_vars.sb.cutoff());
    for (auto const &aa1 : amino_acid::all()) {
      for (auto const &aa2 : amino_acid::all()) {
        cutoff = max(cutoff, params.lj_vars.ss.at({aa1, aa2}).cutoff());
      }
    }
    update.cutoff = cutoff;

#pragma omp single nowait
    st.max_cutoff = max(st.max_cutoff, cutoff);
  }
}

void thread::post_equil_setup() { setup_afm(); }

void thread::setup_afm() {
  if (!params.afm.tips.empty()) {
    auto &vel_eval = eval_vel_afm_forces;
    vel_eval.r = st.r.view();
    vel_eval.t = &st.t;
    vel_eval.afm_force.H1 = params.afm.H1;
    vel_eval.afm_force.H2 = params.afm.H2;
    vel_eval.afm_tips = st.vel_afm_tips.view();
    vel_eval.V = &dyn.V;
    vel_eval.F = dyn.F.view();

    auto &force_eval = eval_force_afm_forces;
    force_eval.afm_tips = st.force_afm_tips.view();
    force_eval.F = dyn.F.view();
  }
}

thread::thread(state &st) : st{st}, params{st.params} {}

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
#pragma omp master
  {
    st.dyn.reset();

    if (params.qa.enabled) {
      qa_finish_processing();
      prepare_nh();
      if (st.ss_spec_crit)
        count_cys_neigh();
    }

    if (st.nl_required) {
      nl_verify();
      if (st.nl_invalid)
        fix_nl();
    }
  };
#pragma omp barrier
}

void thread::fix_nl() {
  switch (params.nl.algorithm) {
  case nl::parameters::CELL:
    nl_cell();
    break;
  case nl::parameters::LEGACY:
    nl_legacy();
    break;
  }

  if (params.pauli.enabled)
    update_pauli_pairs();
  if (params.nat_cont.enabled)
    update_nat_contacts();
  if (params.const_dh.enabled || params.rel_dh.enabled)
    update_dh_pairs();
  if (params.qa.enabled) {
    update_qa_pairs();
    if (st.ss_spec_crit)
      update_cys_neigh();
  }
  if (params.pid.enabled)
    update_pid_bundles();
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
    if (params.pbar.enabled)
      render_pbar();

    if (params.out.enabled)
      make_report();

    if (params.lang.enabled) {
      if (params.lang.type == lang::lang_type::NORMAL)
        lang_step();
      else
        lang_legacy_step();
    }
  }
}
