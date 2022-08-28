#include <cg/simul/thread.h>
#include <omp.h>
namespace cg::simul {

thread_team::thread_team(state &st) {
  this->st = &st;
  num_threads = st.params.gen.num_of_threads;
}

thread &thread_team::fork() {
  thread *thr;
#pragma omp critical
  {
    auto thread_ptr = std::make_unique<thread>(*this);
    thr = threads.emplace_back(std::move(thread_ptr)).get();
  }
  return *thr;
}

thread::thread(thread_team &team) {
  this->team = &team;
  st = team.st;
  params = &st->params;
  tid = omp_get_thread_num();
  eval_slices = set_of_task_slices(team.num_threads);
}

void thread::init_kernels() {
  setup_gen();
  setup_dyn();
  setup_output();
  setup_langevin();
  setup_pbar();
  setup_nl();
  setup_local_rep();
  setup_chir();
  setup_tether();
  setup_angles();
  setup_pauli();
  setup_nat_cont();
  setup_dh();
  setup_qa();
  setup_pid();
  setup_afm();
  setup_walls();
}

void thread::setup_gen() {
  gen = st->gen;
  for (int rep = 0; rep < tid; ++rep)
    gen = gen.spawn();

  if (params->gen.fp_exceptions) {
    feclearexcept(FE_ALL_EXCEPT);
    feenableexcept(FE_INVALID | FE_DIVBYZERO);
  }
}

void thread::setup_dyn() {
  dyn = st->dyn;
#pragma omp critical
  team->forces.push_back(dyn.F);
}

void thread::setup_output() {
  if (st->out_enabled) {
    make_report.struct_every = params->out.struct_every;
    make_report.stats_every = params->out.stats_every;
    make_report.prefix = params->out.prefix;
    make_report.rep = &st->rep;
    make_report.st = st;
    make_report.eval_pid = nullptr;
    make_report.count_pid = std::nullopt;
    make_report.qa = nullptr;
    make_report.nc = nullptr;
    make_report.vel_afm = nullptr;
    make_report.force_afm = nullptr;
  }

  if (params->ckpt.enabled) {
    make_checkpoint.path_fmt = params->ckpt.path_fmt;
    make_checkpoint.every = params->ckpt.every;
    make_checkpoint.st = st;
    make_checkpoint.last_t = &st->ckpt_last_t;
  }

  if (params->gen.dump_data) {
#pragma omp master
    {
      print_raw_data.data_file =
          std::make_shared<std::ofstream>("raw_data.txt");
      print_raw_data.st = st;
    }
  }
}

void thread::setup_langevin() {
  if (st->lang_enabled) {
    auto &step = lang_step;

    step.set_params(params->lang.gamma, st->temperature, params->lang.dt);
    step.gen = &gen;
    step.num_particles = st->num_res;
    step.atype = st->atype;

    step.mass = st->comp_aa_data.mass;
    step.mass_inv = st->mass_inv;
    step.mass_rsqrt = st->mass_rsqrt;

    step.step_idx = &st->step_idx;
    step.t = &st->t;
    step.r = st->r;
    step.v = st->v;
    step.F = st->dyn.F;

    step.y0 = st->y0;
    step.y1 = st->y1;
    step.y2 = st->y2;
    step.y3 = st->y3;
    step.y4 = st->y4;
    step.y5 = st->y5;
    step.true_t = &st->true_t;

    step.noise = st->noise;
  }
}

void thread::setup_pbar() {
  if (st->pbar_enabled) {
    auto &render = render_pbar;
    render.width = params->pbar.width;
    render.total_time = params->gen.total_time;
    render.period_s = params->pbar.update_period.value_in("s");
    render.t = &st->t;
    render.V = &st->dyn.V;
    render.start_clock = &st->pbar_start_clock;
    render.last_clock = &st->pbar_last_clock;
    render.t = &st->t;

    render.is_first = &st->pbar_first_time;
  }
}

void thread::setup_nl() {
  auto &legacy = nl_legacy;
  legacy.pad = params->nl.pad;
  legacy.r = st->r;
  legacy.simul_box = &st->pbc;
  legacy.chain_idx = st->chain_idx;
  legacy.seq_idx = st->seq_idx;
  legacy.num_particles = st->num_res;
  legacy.nl_data = &st->nl;
  legacy.invalid = &st->nl_invalid;
  legacy.cutoff = params->nl.cutoff;
  cutoff = &legacy.cutoff;

  auto &verify = nl_verify;
  verify.r = st->r;
  verify.nl_data = &st->nl;
  verify.simul_box = &st->pbc;
  verify.invalid = &st->nl_invalid;
  verify.first_time = &st->verify_first_time;
  verify.num_particles = st->num_res;
  verify.total_disp = &st->total_disp;
}

void thread::setup_local_rep() {
  if (st->lrep_enabled) {
    auto &eval = eval_lrep_forces;
    eval.F = dyn.F;
    eval.V = &dyn.V;
    eval.depth = params->lrep.depth;
    eval.cutoff = params->gen.repulsive_cutoff;
    eval.r = st->r;
    eval.pairs = st->local_rep_pairs;

    eval_slices.add_static_task(eval);
  }
}

void thread::setup_chir() {
  if (st->chir_enabled) {
    auto &eval = eval_chir_forces;
    eval.e_chi = params->chir.e_chi;
    eval.quads = st->chir_quads;
    eval.V = &dyn.V;
    eval.F = dyn.F;

    eval_slices.add_static_task(eval);
  }
}

void thread::setup_tether() {
  if (params->tether.enabled) {
    auto &eval = eval_tether_forces;
    eval.H1 = params->tether.H1;
    eval.H2 = params->tether.H2;
    eval.def_length = params->tether.def_length;
    eval.r = st->r;
    eval.tethers = st->tether_pairs;
    eval.V = &dyn.V;
    eval.F = dyn.F;

    eval_slices.add_static_task(eval);
  }
}

void thread::setup_angles() {
  if (params->angles.nat_ang.enabled) {
    auto &eval = eval_nat_ang_forces;
    eval.CBA = params->angles.nat_ang.CBA;
    eval.r = st->r;
    eval.angles = st->native_angles;
    eval.V = &dyn.V;
    eval.F = dyn.F;

    eval_slices.add_static_task(eval);
  }

  if (params->angles.heur_ang.enabled) {
    auto &eval = eval_heur_ang_forces;
    eval.r = st->r;
    eval.angles = st->heur_angles;
    eval.V = &dyn.V;
    eval.F = dyn.F;

    for (auto const &heur_pair : aa_heur_pair::all()) {
      for (int d = 0; d <= heur_ang::eval_forces::POLY_DEG; ++d) {
        eval.poly_coeffs[d][(uint8_t)heur_pair] =
            params->angles.heur_ang.coeffs.at(heur_pair).poly[d];
      }
    }

    eval_slices.add_static_task(eval);
  }

  if (params->angles.nat_dih.enabled) {
    auto nat_dih_var = params->angles.nat_dih.variant;
    if (nat_dih_var == "complex") {
      auto &eval = eval_cnd_forces;
      eval.CDA = params->angles.nat_dih.complex.CDA;
      eval.CDB = params->angles.nat_dih.complex.CDB;
      eval.r = st->r;
      eval.dihedrals = st->native_dihedrals;
      eval.V = &dyn.V;
      eval.F = dyn.F;

      eval_slices.add_static_task(eval);
    } else if (nat_dih_var == "simple") {
      auto &eval = eval_snd_forces;
      eval.CDH = params->angles.nat_dih.simple.CDH;
      eval.r = st->r;
      eval.dihedrals = st->native_dihedrals;
      eval.V = &dyn.V;
      eval.F = dyn.F;

      eval_slices.add_static_task(eval);
    }
  }

  if (params->angles.heur_dih.enabled) {
    auto &eval = eval_heur_dih_forces;
    eval.r = st->r;
    eval.V = &dyn.V;
    eval.F = dyn.F;
    eval.dihedrals = st->heur_dihedrals;

    for (auto const &heur_pair : aa_heur_pair::all()) {
      auto idx = (uint8_t)heur_pair;
      auto const &coeffs = params->angles.heur_dih.coeffs.at(heur_pair);
      eval.coeffs.const_[idx] = coeffs.const_;
      eval.coeffs.sin[idx] = coeffs.sin;
      eval.coeffs.cos[idx] = coeffs.cos;
      eval.coeffs.sin2[idx] = coeffs.sin2;
      eval.coeffs.cos2[idx] = coeffs.cos2;
      eval.coeffs.sin_cos[idx] = coeffs.sin_cos;
    }

    eval_slices.add_static_task(eval);
  }
}

void thread::setup_pauli() {
  if (st->pauli_enabled) {
    auto &eval = eval_pauli_forces;
    eval.r_excl = params->gen.repulsive_cutoff;
    eval.depth = params->pauli.depth;
    eval.r = st->r;
    eval.simul_box = &st->pbc;
    eval.pairs = &st->pauli_pairs;
    eval.V = &dyn.V;
    eval.F = dyn.F;

    eval_slices.add_static_task(eval);

    auto &update = update_pauli_pairs;
    update.r = st->r;
    update.simul_box = &st->pbc;
    update.nl = &st->nl;
    update.pairs = &st->pauli_pairs;
    update.r_excl = params->gen.repulsive_cutoff;
  }
}

void thread::setup_nat_cont() {
  if (params->nat_cont.enabled) {
    nl_legacy.all_nat_cont = st->nat_cont_excl;

    auto &eval = eval_nat_cont_forces;
    eval.depth = params->nat_cont.lj_depth;
    eval.simul_box = &st->pbc;
    eval.contacts = &st->cur_native_contacts;
    eval.r = st->r;
    eval.breaking_threshold = params->gen.counting_factor;
    eval.t = &st->t;
    eval.all_contacts = st->all_native_contacts;
    eval.V = &dyn.V;
    eval.F = dyn.F;
    eval.cutoff = params->nl.cutoff;
    eval.num_changed = &st->num_changed;

    auto get_disul_force = [&](force_spec spec) -> disulfide_force {
      disulfide_force res;
      if (spec.variant == "harmonic")
        res.force = (harmonic)spec.harmonic.value();
      else if (spec.variant == "lj")
        res.force = (lj)spec.lj.value();
      else
        throw std::runtime_error("invalid force spec for disulfide force");
      return res;
    };

    if (params->nat_cont.ss_force.has_value())
      eval.disulfide = get_disul_force(params->nat_cont.ss_force.value());
    else
      eval.disulfide = std::nullopt;

    eval_slices.add_static_task(eval);

    auto &update = update_nat_contacts;
    update.r = st->r;
    update.simul_box = &st->pbc;
    update.nl = &st->nl;
    update.all_contacts = st->all_native_contacts;
    update.contacts = &st->cur_native_contacts;
    update.cutoff = params->nl.cutoff;

    if (params->out.enabled)
      make_report.nc = &eval;
  }
}

void thread::setup_dh() {
  auto &update = update_dh_pairs;
  update.r = st->r;
  update.simul_box = &st->pbc;
  update.nl = &st->nl;
  update.pairs = &st->dh_pairs;
  update.atype = st->atype;
  update.cutoff = params->dh.cutoff;
  for (auto const &aa : amino_acid::all()) {
    update.q[(uint8_t)aa] = st->comp_aa_data.charge[(uint8_t)aa];
  }

  *cutoff = max(*cutoff, update.cutoff);

  auto dh_var = params->dh.variant;
  if (dh_var == "constant") {
    auto &eval = eval_const_dh_forces;
    eval.set_V_factor(params->dh.const_dh.permittivity);
    eval.screen_dist_inv = 1.0 / params->dh.screening_dist;
    eval.r = st->r;
    eval.simul_box = &st->pbc;
    eval.es_pairs = st->dh_pairs;
    eval.V = &dyn.V;
    eval.F = dyn.F;
    eval.cutoff = params->dh.cutoff;

    eval_slices.add_static_task(eval);
  } else if (dh_var == "relative") {
    auto &eval = eval_rel_dh_forces;
    eval.set_V_factor(params->dh.rel_dh.perm_factor);
    eval.screen_dist_inv = 1.0 / params->dh.screening_dist;
    eval.r = st->r;
    eval.simul_box = &st->pbc;
    eval.es_pairs = st->dh_pairs;
    eval.V = &dyn.V;
    eval.F = dyn.F;
    eval.cutoff = params->dh.cutoff;

    eval_slices.add_static_task(eval);
  }
}

void thread::setup_qa() {
  if (params->qa.enabled) {
    auto &prep_nh = prepare_nh;
    prep_nh.r = st->r;
    prep_nh.num_particles = st->num_res;
    prep_nh.n = st->n;
    prep_nh.h = st->h;
    prep_nh.prev = st->prev;
    prep_nh.next = st->next;
    prep_nh.simul_box = &st->pbc;

    auto &loop = qa_loop_over_candidates;
    loop.r = st->r;
    loop.n = st->n;
    loop.h = st->h;
    loop.simul_box = &st->pbc;
    loop.atype = st->atype;
    loop.sync = st->sync_values;
    loop.free_pairs = &st->qa_free_pairs;
    loop.candidates = &st->qa_candidates;
    loop.total_disp = &st->total_disp;

    loop.min_abs_cos_hr = params->qa.min_cos_hr;
    loop.min_abs_cos_hh = params->qa.min_cos_hh;
    loop.max_cos_nr = params->qa.max_cos_nr;
    loop.formation_tolerance = params->qa.formation_tolerance;
    loop.disulfide_special_criteria = st->ss_spec_crit;

    loop.rep_cutoff = params->gen.repulsive_cutoff;
    loop.rep_depth = params->pauli.depth;
    loop.F = dyn.F;
    loop.V = &dyn.V;

    for (auto const &aa : amino_acid::all())
      loop.ptype[(uint8_t)aa] = st->comp_aa_data.ptype[(uint8_t)aa];

    eval_slices.add_static_task(loop);

    auto &fin_proc = qa_finish_processing;
    fin_proc.candidates = &st->qa_candidates;
    fin_proc.sync = st->sync_values;
    fin_proc.t = &st->t;
    fin_proc.contacts = &st->qa_contacts;
    fin_proc.free_pairs = &st->qa_free_pairs;
    fin_proc.removed = &st->qa_removed;
    fin_proc.num_contacts = &st->num_qa_contacts;
    fin_proc.disulfide_special_criteria = st->ss_spec_crit;

    auto &proc_cont = process_qa_contacts;
    proc_cont.saturation_speed = (real)1.0 / params->qa.phase_dur;
    proc_cont.dt = params->lang.dt;
    proc_cont.set_factor(params->gen.counting_factor);
    proc_cont.t = &st->t;
    proc_cont.sync = st->sync_values;
    proc_cont.contacts = &st->qa_contacts;
    proc_cont.simul_box = &st->pbc;
    proc_cont.r = st->r;
    proc_cont.free_pairs = &st->qa_free_pairs;
    proc_cont.removed = &st->qa_removed;
    proc_cont.V = &dyn.V;
    proc_cont.F = dyn.F;
    proc_cont.disulfide_special_criteria = st->ss_spec_crit;
    proc_cont.fixed_cutoff = params->nl.cutoff;

    auto get_force = [&](force_spec spec) -> sink_lj {
      if (spec.variant == "lj") {
        return sink_lj((lj)spec.lj.value());
      } else if (spec.variant == "sink lj") {
        auto sink_lj_ = spec.sink_lj.value_or(sink_lj_specs());
        auto lj_ = spec.lj.value();

        if (!sink_lj_.depth.has_value())
          sink_lj_.depth = lj_.depth;
        if (!sink_lj_.r_low.has_value())
          sink_lj_.r_low = params->gen.repulsive_cutoff;
        if (!sink_lj_.r_high.has_value())
          sink_lj_.r_high = lj_.r_min;

        return sink_lj_;
      } else {
        throw std::runtime_error("invalid force spec for QA potential");
      }
    };

    auto get_ss_force = [&](ss_force_spec spec, amino_acid aa1,
                            amino_acid aa2) -> sink_lj {
      if (spec.variant == "lj") {
        return sink_lj((lj)spec.lj->ss_specs.at({aa1, aa2}));
      } else if (spec.variant == "sink lj") {
        auto sink_lj_ = spec.sink_lj.has_value()
                            ? spec.sink_lj->ss_specs.at({aa1, aa2})
                            : sink_lj_specs();
        auto lj_ = spec.lj->ss_specs.at({aa1, aa2});

        if (!sink_lj_.depth.has_value())
          sink_lj_.depth = lj_.depth;
        if (!sink_lj_.r_low.has_value())
          sink_lj_.r_low = params->gen.repulsive_cutoff;
        if (!sink_lj_.r_high.has_value())
          sink_lj_.r_high = lj_.r_min;

        return sink_lj_;
      } else {
        throw std::runtime_error("invalid force spec for QA potential");
      }
    };

    proc_cont.ljs = vect::vector<sink_lj>(qa::contact_type::NUM_TYPES);
    proc_cont.ljs[(int16_t)qa::contact_type::BACK_BACK()] =
        get_force(params->qa.bb);
    proc_cont.ljs[(int16_t)qa::contact_type::BACK_SIDE()] =
        proc_cont.ljs[(int16_t)qa::contact_type::SIDE_BACK()] =
            get_force(params->qa.bs);

    for (auto const &aa1 : amino_acid::all()) {
      for (auto const &aa2 : amino_acid::all()) {
        proc_cont.ljs[(int16_t)qa::contact_type::SIDE_SIDE(aa1, aa2)] =
            get_ss_force(params->qa.ss, aa1, aa2);
      }
    }

    for (auto const &ctype : qa::contact_type::all())
      loop.req_min_dist[(int16_t)ctype] =
          proc_cont.ljs[(int16_t)ctype].r_high();

    eval_slices.add_static_task(proc_cont);

    auto &update = update_qa_pairs;
    update.r = st->r;
    update.simul_box = &st->pbc;
    update.nl = &st->nl;
    update.pairs = &st->qa_free_pairs;
    update.chain_idx = st->chain_idx;
    update.seq_idx = st->seq_idx;
    update.include4 = params->qa.include4;

    if (params->qa.disulfide.has_value()) {
      if (st->ss_spec_crit) {
        auto const &spec_crit_params = params->qa.disulfide->spec_crit;

        loop.part_of_ssbond = st->part_of_ssbond;
        loop.disulfide_special_criteria = st->ss_spec_crit;
        loop.neigh = st->neigh_count;

        fin_proc.part_of_ssbond = st->part_of_ssbond;

        auto get_disul_force = [&](force_spec spec) -> disulfide_force {
          disulfide_force res;
          if (spec.variant == "harmonic")
            res.force = (harmonic)spec.harmonic.value();
          else if (spec.variant == "lj")
            res.force = (lj)spec.lj.value();
          else
            throw std::runtime_error("invalid force spec for disulfide force");
          return res;
        };

        proc_cont.disulfide = get_disul_force(params->qa.disulfide->force);
        proc_cont.ss_def_dist = params->qa.disulfide->spec_crit.def_dist;
        proc_cont.ss_dist_max_div =
            params->qa.disulfide->spec_crit.max_dist_dev;
        proc_cont.neigh = st->neigh_count;
        proc_cont.max_neigh_count = spec_crit_params.max_neigh_count;
        proc_cont.part_of_ssbond = st->part_of_ssbond;

        auto &update_cys = update_cys_neigh;
        update_cys.neigh = &st->qa_cys_neigh;
        update_cys.neigh_radius = spec_crit_params.neigh_radius;
        update_cys.nl = &st->nl;
        update_cys.simul_box = &st->pbc;
        update_cys.r = st->r;

        auto &count_cys = count_cys_neigh;
        count_cys.simul_box = &st->pbc;
        count_cys.neigh_radius = spec_crit_params.neigh_radius;
        count_cys.neigh = &st->qa_cys_neigh;
        count_cys.neigh_count = st->neigh_count;
        count_cys.cys_indices = st->cys_indices;
        count_cys.r = st->r;

        auto &reset_cys = reset_cys_neigh;
        reset_cys.cys_indices = count_cys.cys_indices;
        reset_cys.neigh_count = count_cys.neigh_count;
      }
    }

    real max_req_dist = 0.0;
    for (auto const &x : loop.req_min_dist)
      max_req_dist = max(max_req_dist, x);

    update.max_formation_min_dist = loop.max_req_dist = max_req_dist;

    update.cutoff = params->nl.cutoff;

    if (params->out.enabled)
      make_report.qa = &proc_cont;
  }
}

void thread::setup_pid() {
  if (params->pid.enabled) {
    auto &eval = eval_pid_forces;
    eval.V = &dyn.V;
    eval.F = dyn.F;
    eval.total_disp = &st->total_disp;
    eval.active_thr = params->gen.counting_factor;
    eval.atype = st->atype;

    eval.bb_plus_lam.version() = params->pid.lambda_variant;
    eval.bb_plus_lam.alpha() = params->pid.bb_plus_lambda.alpha;
    eval.bb_plus_lam.psi_0() = params->pid.bb_plus_lambda.psi_0;

    eval.bb_minus_lam.version() = params->pid.lambda_variant;
    eval.bb_minus_lam.alpha() = params->pid.bb_minus_lambda.alpha;
    eval.bb_minus_lam.psi_0() = params->pid.bb_minus_lambda.psi_0;

    eval.ss_lam.version() = params->pid.lambda_variant;
    eval.ss_lam.alpha() = params->pid.ss_lambda.alpha;
    eval.ss_lam.psi_0() = params->pid.ss_lambda.psi_0;

    auto get_force = [&](force_spec spec) -> sink_lj {
      if (spec.variant == "lj") {
        return sink_lj((lj)spec.lj.value());
      } else if (spec.variant == "sink lj") {
        auto sink_lj_ = spec.sink_lj.value_or(sink_lj_specs());
        auto lj_ = spec.lj.value();

        if (!sink_lj_.depth.has_value())
          sink_lj_.depth = lj_.depth;
        if (!sink_lj_.r_low.has_value())
          sink_lj_.r_low = 0;
        if (!sink_lj_.r_high.has_value())
          sink_lj_.r_high = lj_.r_min;

        return sink_lj_;
      } else {
        throw std::runtime_error("invalid force spec for PID potential");
      }
    };

    auto get_ss_force = [&](ss_force_spec spec, amino_acid aa1,
                            amino_acid aa2) -> sink_lj {
      if (spec.variant == "lj") {
        return sink_lj((lj)spec.lj->ss_specs.at({aa1, aa2}));
      } else if (spec.variant == "sink lj") {
        auto sink_lj_ = spec.sink_lj.has_value()
                            ? spec.sink_lj->ss_specs.at({aa1, aa2})
                            : sink_lj_specs();
        auto lj_ = spec.lj->ss_specs.at({aa1, aa2});

        if (!sink_lj_.depth.has_value())
          sink_lj_.depth = lj_.depth;
        if (!sink_lj_.r_low.has_value())
          sink_lj_.r_low = 0;
        if (!sink_lj_.r_high.has_value())
          sink_lj_.r_high = lj_.r_min;

        return sink_lj_;
      } else {
        throw std::runtime_error("invalid force spec for PID potential");
      }
    };

    eval.bb_plus_lj = get_force(params->pid.bb_plus_force);
    eval.bb_minus_lj = get_force(params->pid.bb_minus_force);

    eval.ss_ljs =
        vect::vector<sink_lj>(amino_acid::NUM_TYPES * amino_acid::NUM_TYPES);
    for (auto const &aa1 : amino_acid::all()) {
      auto idx1 = (uint16_t)(uint8_t)aa1;
      for (auto const &aa2 : amino_acid::all()) {
        auto idx2 = (uint16_t)(uint8_t)aa2;
        auto ss_idx = idx1 * (uint16_t)amino_acid::NUM_TYPES + idx2;
        eval.ss_ljs[ss_idx] = get_ss_force(params->pid.ss_force, aa1, aa2);
      }
    }

    eval.r = st->r;
    eval.simul_box = &st->pbc;
    eval.bundles = &st->pid_bundles;
    eval.prev = st->prev;
    eval.next = st->next;
    eval.V = &dyn.V;
    eval.F = dyn.F;

    eval_slices.add_static_task(eval);

    auto &update = update_pid_bundles;
    update.r = st->r;
    update.prev = st->prev;
    update.next = st->next;
    update.atype = st->atype;
    update.simul_box = &st->pbc;
    update.nl = &st->nl;
    update.bundles = &st->pid_bundles;
    update.chain_idx = st->chain_idx;
    update.seq_idx = st->seq_idx;
    update.include4 = params->pid.include4;

    eval.cutoff = update.cutoff = params->nl.cutoff;

    if (params->out.enabled) {
      make_report.eval_pid = &eval;
      auto &count = make_report.count_pid.emplace();
      count.eval = &eval;
      count.chain_idx = st->chain_idx;
      count.count_factor = params->gen.counting_factor;
    }
  }
}

void thread::setup_afm() {
  if (st->afm_enabled) {
    auto &eval_vafm = eval_vel_afm_forces;
    eval_vafm.r = st->r;
    eval_vafm.F = dyn.F;
    eval_vafm.V = &dyn.V;
    eval_vafm.t = &st->t;
    eval_vafm.afm_force.H1 = params->afm.tip_params.vel_afm.H1;
    eval_vafm.afm_force.H2 = params->afm.tip_params.vel_afm.H2;
    eval_vafm.afm_force.nat_r = 0;
    eval_vafm.afm_tips = st->vel_afm_tips;

    eval_slices.add_static_task(eval_vafm);

    auto &eval_afm = eval_force_afm_forces;
    eval_afm.F = dyn.F;
    eval_afm.afm_tips = st->force_afm_tips;

    eval_slices.add_static_task(eval_afm);

    if (params->out.enabled) {
      make_report.vel_afm = &eval_vafm;
      make_report.force_afm = &eval_afm;
    }
  }
}

void thread::setup_solid_walls() {
  if (st->solid_walls_enabled) {
    auto &eval = eval_solid_wall_forces;
    eval.min_dist = params->sbox.walls.threshold;
    eval.walls = st->solid_walls;
    eval.wall_F = dyn.solid_wall_F;
    eval.depth = params->sbox.walls.solid_wall.depth;
    eval.r = st->r;
    eval.F = dyn.F;

    eval_slices.add_static_task(eval);
  }
}

void thread::setup_harmonic_walls() {
  if (st->harmonic_walls_enabled) {
    auto &free = hw_eval_free;
    free.depth = params->sbox.walls.solid_wall.depth;
    free.walls = st->harmonic_walls;
    free.min_dist = params->sbox.walls.threshold;
    free.is_connected = st->is_connected_to_wall;
    free.wall_F = dyn.harmonic_wall_F;
    free.r = st->r;
    free.F = dyn.F;
    free.V = &dyn.V;

    eval_slices.add_static_task(free);

    auto &conn = hw_eval_conn;
    conn.walls = st->harmonic_walls;
    conn.wall_F = dyn.harmonic_wall_F;
    conn.conns = st->harmonic_conns;
    conn.r = st->r;
    conn.HH1 = params->sbox.walls.harmonic_wall.HH1;
    conn.F = dyn.F;
    conn.V = &dyn.V;

    eval_slices.add_static_task(conn);
  }
}

void thread::setup_lj_walls() {
  if (st->lj_walls_enabled) {
    auto &sift = ljw_sift_free;
    sift.r = st->r;
    sift.walls = st->lj_walls;
    sift.is_connected = st->is_connected_to_wall;
    sift.min_dist = params->sbox.walls.threshold;
    sift.removed = &st->ljw_removed;
    sift.candidates = &st->ljw_candidates;
    sift.wall_F = dyn.lj_wall_F;
    sift.force.depth = params->sbox.walls.lj_wall.depth;
    sift.force.r_min = sift.min_dist;
    sift.F = dyn.F;
    sift.V = &dyn.V;

    eval_slices.add_static_task(sift);

    auto &eval = ljw_eval_conn;
    eval.min_dist = sift.min_dist;
    eval.force.depth() = sift.force.depth;
    eval.force.r_min() = sift.force.r_min;
    eval.walls = st->lj_walls;
    eval.removed = sift.removed;
    eval.conns = &st->ljw_conns;
    eval.r = sift.r;
    eval.is_connected = st->is_connected_to_wall;
    eval.breaking_dist = params->sbox.walls.lj_wall.breaking_dist_factor *
                         pow(0.5, 1.0 / 6.0) * eval.min_dist;
    eval.saturation_diff =
        params->lang.dt / params->sbox.walls.lj_wall.cycle_dur;
    eval.wall_F = sift.wall_F;
    eval.F = sift.F;
    eval.V = sift.V;

    eval_slices.add_static_task(eval);

    auto &proc = ljw_proc_cand;
    proc.walls = eval.walls;
    proc.is_connected = eval.is_connected;
    proc.conns = eval.conns;
    proc.removed = eval.removed;
    proc.candidates = sift.candidates;
    proc.r = eval.r;
  }
}

void thread::setup_walls() {
  setup_solid_walls();
  setup_harmonic_walls();
  setup_lj_walls();
  dyn = st->dyn;

  log_wall_forces_enabled = st->walls[POS_Z] && st->walls[NEG_Z];
  if (log_wall_forces_enabled) {
    auto &log = log_wall_forces;
    log.walls = st->walls;
    log.neg_z = st->walls[NEG_Z];
    log.pos_z = st->walls[POS_Z];
    log.t = &st->t;
    log.avg_z_force = &st->avg_z_force;
  }
}
} // namespace cg::simul