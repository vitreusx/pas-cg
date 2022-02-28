#include "simul/simulation.h"
#include <iostream>
using namespace cg::simul;

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <fenv.h>

void simulation::print_help(char **argv) {
  std::cout << "Usage: " << argv[0] << " [--help|file]" << '\n';
  std::cout << '\t' << "file" << '\t' << "parameter file" << '\n';
  std::cout << '\t' << "--help" << '\t' << "display this help message";
}

void simulation::parse_args(int argc, char **argv) {
  std::vector<std::string> argv_s(argc);
  for (int idx = 0; idx < argc; ++idx)
    argv_s[idx] = argv[idx];

  if (argc == 2 && argv_s[1] == "--help") {
    print_help(argv);
    exit(EXIT_SUCCESS);
  } else if (argc > 2) {
    print_help(argv);
    exit(EXIT_FAILURE);
  } else {
    if (argc > 1)
      param_path = argv_s[1];
  }
}

void simulation::load_parameters() {
  using namespace ioxx::xyaml;
  auto params_yml = node::import("data/default/inputfile.yml");
  //  if (param_path) {
  //    auto overrides = ioxx::xyaml_node::from_path(param_path.value());
  //    params_yml = merge(params, overrides);
  //  }

  params_yml >> params;
}

void simulation::general_setup() {
  omp_set_num_threads((int)params.gen.num_of_threads);
  gen = rand_gen(params.gen.seed);

  if (params.gen.debug_mode) {
    feclearexcept(FE_ALL_EXCEPT);
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
  }
}

void simulation::load_model() {
  auto &model_file_v = params.input.source;
  if (std::holds_alternative<pdb_file>(model_file_v)) {
    auto const &file = std::get<pdb_file>(model_file_v);
    model = file.to_model();
  } else {
    auto &file = std::get<seq_file>(model_file_v);
    model = std::move(file.model);
  }

  if (params.input.morph_into_saw.has_value()) {
    auto &saw_p = params.input.morph_into_saw.value();
    if (saw_p.perform)
      model.morph_into_saw(gen, saw_p.bond_distance, saw_p.residue_density,
                           saw_p.infer_box);
  }
}

void simulation::compile_model() {
  int res_idx = 0;
  for (auto const &res : model.residues)
    res_map[res.get()] = res_idx++;

  num_res = (int)model.residues.size();

  r = nitro::vector<vec3r>(num_res);
  for (auto const &res : model.residues)
    r[res_map.at(res.get())] = res->pos;

  atype = nitro::vector<amino_acid>(num_res);
  for (auto const &res : model.residues)
    atype[res_map.at(res.get())] = res->type;

  comp_aa_data = compiled_aa_data(params.aa_data);

  box.cell = model.model_box.cell;
  box.cell_inv = model.model_box.cell_inv;

  prev = next = chain_idx = seq_idx = nitro::vector<int>(num_res, -1);
  for (auto const &res : model.residues) {
    if (res->seq_idx > 0) {
      auto prev_ptr = res->parent->residues[res->seq_idx - 1];
      prev[res_map.at(res.get())] = res_map.at(prev_ptr);
    }

    auto num_res_in_chain = (int)res->parent->residues.size();
    if (res->seq_idx < num_res_in_chain - 1) {
      auto next_ptr = res->parent->residues[res->seq_idx + 1];
      next[res_map.at(res.get())] = res_map.at(next_ptr);
    }

    chain_idx[res_map.at(res.get())] = res->parent->chain_idx;
    seq_idx[res_map.at(res.get())] = res->seq_idx;
  }
}

void simulation::setup_dyn() {
  t = 0;
  V = 0;
  dyn = dynamics(num_res);
}

void simulation::setup_output() {
  if (params.out.enabled) {
    auto &make_rep = ker.make_report;
    report.first_time = true;
    report.ord = 0;
    report.output_dir = params.out.output_dir;
    make_rep.period = params.out.period;
    make_rep.t = &t;
    report_last_t = std::numeric_limits<real>::lowest();
    make_rep.last_t = &report_last_t;
    make_rep.hooks = &hooks;
    make_rep.state = &report;

    auto &export_pdb = ker.export_pdb;
    export_pdb.ref_model = &model;
    export_pdb.r = r.view();
    export_pdb.res_map = &res_map;
    hooks.emplace_back(&export_pdb);

    auto &add_stats = ker.add_stats;
    add_stats.t = &t;
    add_stats.V = &dyn.V;
    add_stats.v = v.view();
    hooks.emplace_back(&add_stats);

    auto &add_structure = ker.add_structure;
    add_structure.res_map = &res_map;
    add_structure.r = r.view();
    add_structure.model = &model;
    hooks.emplace_back(&add_structure);

    auto& report_gyr = ker.report_gyr;
    report_gyr.r = r.view();
    hooks.emplace_back(&report_gyr);
  }
}

void simulation::setup_langevin() {
  if (params.lang.enabled) {
    auto &step = ker.lang_step;
    step.temperature = params.lang.temperature;
    step.t = &t;
    step.dt = params.lang.dt;
    step.gamma_factor = params.lang.gamma;
    step.gen = &gen;
    step.mass = comp_aa_data.mass.view();
    step.num_particles = num_res;
    step.atype = atype.view();
    step.r = r.view();

    mass_inv = nitro::vector<real>(step.mass.size());
    for (int aa_idx = 0; aa_idx < mass_inv.size(); ++aa_idx)
      mass_inv[aa_idx] = 1.0 / step.mass[aa_idx];
    step.mass_inv = mass_inv.view();

    mass_rsqrt = nitro::vector<real>(step.mass.size());
    for (int aa_idx = 0; aa_idx < mass_rsqrt.size(); ++aa_idx)
      mass_rsqrt[aa_idx] = sqrt(mass_inv[aa_idx]);
    step.mass_rsqrt = mass_rsqrt.view();

    v = nitro::vector<vec3r>(num_res, vec3r::Zero());
    step.v = v.view();

    y0 = y1 = y2 = y3 = y4 = y5 = nitro::vector<vec3sr>(num_res);
    for (int idx = 0; idx < num_res; ++idx)
      y0[idx] = r[idx];
    step.y0 = y0.view();
    step.y1 = y1.view();
    step.y2 = y2.view();
    step.y3 = y3.view();
    step.y4 = y4.view();
    step.y5 = y5.view();

    true_t = t;
    step.true_t = &true_t;

    step.F = dyn.F.view();
  }
}

void simulation::setup_pbar() {
  if (params.pbar.enabled) {
    auto &render = ker.render_pbar;
    render.width = params.pbar.width;
    render.total_time = params.gen.total_time;
    render.period = params.pbar.update_period;
    render.t = &t;
    render.V = &dyn.V;
    render.start_wall_time = &start_wall_time;
    render.t = &t;
    render.last_t = &pbar_last_t;

    pbar_first_time = true;
    render.is_first = &pbar_first_time;
  }
}

void simulation::setup_nl() {
  nl_invalid = true;
  nl_required = false;
  max_cutoff = 0.0;

  if (params.nl.algorithm == nl::parameters::LEGACY) {
    auto &legacy = ker.nl_legacy;
    legacy.pad_factor = params.nl.pad_factor;
    legacy.r = r.view();
    legacy.box = &box;
    legacy.t = &t;
    legacy.chain_idx = chain_idx.view();
    legacy.seq_idx = seq_idx.view();
    legacy.num_particles = num_res;
    legacy.data = &nl;
    legacy.data->orig_r = nitro::vector<vec3r>(num_res);
    legacy.invalid = &nl_invalid;
    legacy.max_cutoff = &max_cutoff;
  } else {
    auto &cell = ker.nl_cell;
    cell.pad_factor = params.nl.pad_factor;
    cell.r = r.view();
    cell.box = &box;
    cell.t = &t;
    cell.chain_idx = chain_idx.view();
    cell.seq_idx = seq_idx.view();
    cell.num_particles = num_res;
    cell.data = &nl;
    cell.data->orig_r = nitro::vector<vec3r>(num_res);
    cell.invalid = &nl_invalid;
    cell.max_cutoff = &max_cutoff;
    res_cell_idx = reordered_idx = nitro::vector<int>(num_res);
    cell.res_cell_idx = res_cell_idx.view();
    cell.reordered_idx = reordered_idx.view();
    cell.num_res_in_cell = &num_res_in_cell;
    cell.cell_offset = &cell_offset;
    cell.all_pairs = &all_pairs;
  }

  auto &verify = ker.nl_verify;
  verify.r = r.view();
  verify.data = &nl;
  verify.box = &box;
  verify.invalid = &nl_invalid;
  verify_first_time = true;
  verify.first_time = &verify_first_time;
  verify.num_particles = num_res;
}

void simulation::setup_chir() {
  if (params.chir.enabled) {
    auto &eval = ker.eval_chir_forces;
    eval.e_chi = params.chir.e_chi;

    for (auto const &dihedral : model.dihedrals) {
      auto i1 = res_map[dihedral.res1], i2 = res_map[dihedral.res2],
           i3 = res_map[dihedral.res3], i4 = res_map[dihedral.res4];

      vec3r nat_r1 = dihedral.res1->pos, nat_r2 = dihedral.res2->pos,
            nat_r3 = dihedral.res3->pos, nat_r4 = dihedral.res4->pos;

      auto nat_r12 = nat_r2 - nat_r1, nat_r23 = nat_r3 - nat_r2,
           nat_r34 = nat_r4 - nat_r3;

      auto nat_factor = ipow<3>(norm_inv(nat_r23));
      auto nat_chir = dot(nat_r12, cross(nat_r23, nat_r34)) * nat_factor;

      chir_quads.emplace_back(i1, i2, i3, i4, nat_factor, nat_chir);
    }
    eval.quads = chir_quads.view();
  }
}

void simulation::setup_tether() {
  if (params.tether.enabled) {
    auto &eval = ker.eval_tether_forces;
    eval.H1 = params.tether.H1;
    eval.H2 = params.tether.H2;
    eval.def_length = params.tether.def_length;
    eval.r = r.view();

    for (auto const &tether : model.tethers) {
      auto i1 = res_map[tether.res1], i2 = res_map[tether.res2];
      auto nat_dist = (real)tether.length.value_or(eval.def_length);
      tether_pairs.emplace_back(i1, i2, nat_dist);
    }
    eval.tethers = tether_pairs.view();
  }
}

void simulation::setup_nat_ang() {
  if (params.nat_ang.enabled) {
    auto &eval = ker.eval_nat_ang_forces;
    eval.k = params.nat_ang.k;
    eval.r = r.view();

    for (auto const &angle : model.angles) {
      if (angle.theta.has_value()) {
        auto i1 = res_map[angle.res1], i2 = res_map[angle.res2],
             i3 = res_map[angle.res3];
        auto nat_theta = (real)angle.theta.value();

        native_angles.emplace_back(i1, i2, i3, nat_theta);
      }
    }
    eval.angles = native_angles.view();
  }
}

void simulation::setup_nat_dih() {
  for (auto const &dihedral : model.dihedrals) {
    if (dihedral.phi.has_value()) {
      auto i1 = res_map[dihedral.res1], i2 = res_map[dihedral.res2],
           i3 = res_map[dihedral.res3], i4 = res_map[dihedral.res4];
      auto nat_phi = (real)dihedral.phi.value();

      native_dihedrals.emplace_back(i1, i2, i3, i4, nat_phi);
    }
  }

  if (params.cnd.enabled) {
    auto &eval = ker.eval_cnd_forces;
    eval.CDA = params.cnd.CDA;
    eval.CDB = params.cnd.CDB;
    eval.r = r.view();
    eval.dihedrals = native_dihedrals.view();
  }

  if (params.snd.enabled) {
    auto &eval = ker.eval_snd_forces;
    eval.CDH = params.snd.CDH;
    eval.r = r.view();
    eval.dihedrals = native_dihedrals.view();
  }
}

void simulation::setup_pauli() {
  if (params.pauli.enabled) {
    auto &eval = ker.eval_pauli_forces;
    eval.r_excl = params.pauli.r_excl;
    eval.depth = params.pauli.depth;
    eval.r = r.view();
    eval.box = &box;
    eval.pairs = &pauli_pairs;

    auto &update = ker.update_pauli_pairs;
    update.r = r.view();
    update.box = &box;
    update.nl = &nl;
    update.pairs = &pauli_pairs;

    nl_required = true;
    max_cutoff = max(max_cutoff, (real)params.pauli.r_excl);
  }
}

void simulation::setup_nat_cont() {
  if (params.nat_cont.enabled) {
    for (auto const &cont : model.contacts) {
      if (cont.type != input::model::NAT_SS) {
        auto i1 = res_map[cont.res1], i2 = res_map[cont.res2];
        auto nat_dist = (real)cont.length;
        all_native_contacts.emplace_back(i1, i2, nat_dist);
        native_contact_exclusions.emplace_back(i1, i2, 0.0);
      }
    }

    ker.nl_legacy.all_nat_cont = native_contact_exclusions.view();
    ker.nl_cell.all_nat_cont = native_contact_exclusions.view();

    auto &eval = ker.eval_nat_cont_forces;
    eval.depth = params.nat_cont.lj_depth;
    eval.box = &box;
    eval.contacts = &cur_native_contacts;
    eval.r = r.view();

    auto &update_ref = ker.update_nat_contacts;
    update_ref.r = r.view();
    update_ref.box = &box;
    update_ref.nl = &nl;
    update_ref.all_contacts = all_native_contacts.view();
    update_ref.contacts = &cur_native_contacts;

    real cutoff = 0.0;
    for (int cont_idx = 0; cont_idx < all_native_contacts.size(); ++cont_idx) {
      auto cont_dist = all_native_contacts[cont_idx].nat_dist();
      cutoff = max(cutoff, lj::compute_cutoff(cont_dist));
    }
    nl_required = true;
    max_cutoff = max(max_cutoff, cutoff);

    if (params.out.enabled) {
      auto& report_nc = ker.report_nc_stuff;
      report_nc.all_nat_conts = all_native_contacts.view();
      report_nc.params = &params.nat_cont;
      report_nc.r = r.view();
      hooks.push_back(&report_nc);
    }
  }
}

void simulation::setup_dh() {
  auto &update = ker.update_dh_pairs;
  update.r = r.view();
  update.box = &box;
  update.nl = &nl;
  update.pairs = &dh_pairs;
  update.atype = atype.view();
  for (auto const &aa : amino_acid::all())
    update.q[(uint8_t)aa] = comp_aa_data.charge[(uint8_t)aa];

  if (params.const_dh.enabled) {
    auto &eval = ker.eval_const_dh_forces;
    eval.set_V_factor(params.const_dh.permittivity);
    eval.screen_dist_inv = 1.0 / params.const_dh.screening_dist;
    eval.r = r.view();
    eval.box = &box;
    eval.es_pairs = &dh_pairs;
    update.cutoff = 2.0 * params.const_dh.screening_dist;
  }

  if (params.rel_dh.enabled) {
    auto &eval = ker.eval_rel_dh_forces;
    eval.set_V_factor(params.rel_dh.perm_factor);
    eval.screen_dist_inv = 1.0 / params.rel_dh.screening_dist;
    eval.r = r.view();
    eval.box = &box;
    eval.es_pairs = &dh_pairs;
    update.cutoff = 2.0 * params.rel_dh.screening_dist;
  }

  if (params.const_dh.enabled || params.rel_dh.enabled) {
    nl_required = true;
    max_cutoff = max(max_cutoff, update.cutoff);
  }
}

void simulation::setup_qa() {
  if (params.qa.enabled) {
    sync_values = nitro::vector<sync_data>(num_res);
    for (int res_idx = 0; res_idx < num_res; ++res_idx) {
      auto lim = comp_aa_data.base_sync[(uint8_t)atype[res_idx]];
      sync_values[res_idx] = sync_data(
          lim.back(), lim.side_all(), lim.side_polar(), lim.side_hydrophobic());
    }

    n = h = nitro::vector<vec3r>(num_res);
    auto &prep_nh = ker.prepare_nh;
    prep_nh.r = r.view();
    prep_nh.num_particles = num_res;
    prep_nh.n = n.view();
    prep_nh.h = h.view();
    prep_nh.prev = prev.view();
    prep_nh.next = next.view();

    auto &sift = ker.sift_qa_candidates;
    sift.r = r.view();
    sift.atype = atype.view();
    sift.box = &box;
    sift.free_pairs = &qa_free_pairs;
    sift.sync = sync_values.view();
    sift.candidates = &qa_candidates;
    sift.n = n.view();
    sift.h = h.view();

    sift.min_abs_cos_hr = params.qa.min_cos_hr;
    sift.min_abs_cos_hh = params.qa.min_cos_hh;
    sift.max_cos_nr = params.qa.max_cos_nr;

    sift.req_min_dist[(int16_t)qa::contact_type::BACK_BACK()] =
        params.lj_variants.bb.r_min();
    sift.req_min_dist[(int16_t)qa::contact_type::BACK_SIDE()] =
        params.lj_variants.bs.r_min();
    sift.req_min_dist[(int16_t)qa::contact_type::SIDE_BACK()] =
        params.lj_variants.sb.r_min();

    for (auto const &aa1 : amino_acid::all()) {
      for (auto const &aa2 : amino_acid::all()) {
        sift.req_min_dist[(int16_t)qa::contact_type::SIDE_SIDE(aa1, aa2)] =
            params.lj_variants.ss[{aa1, aa2}].r_max();
      }
    }

    for (auto const &aa : amino_acid::all())
      sift.ptype[(uint8_t)aa] = comp_aa_data.ptype[(uint8_t)aa];

    auto &proc_cand = ker.process_qa_candidates;
    proc_cand.candidates = &qa_candidates;
    proc_cand.contacts = &qa_contacts;
    proc_cand.sync = sync_values.view();
    proc_cand.free_pairs = &qa_free_pairs;
    proc_cand.t = &t;

    auto &proc_cont = ker.process_qa_contacts;
    proc_cont.cycle_time = params.qa.phase_dur;
    proc_cont.cycle_time_inv = 1.0 / proc_cont.cycle_time;
    proc_cont.set_factor(params.qa.breaking_factor);
    proc_cont.t = &t;
    proc_cont.sync = sync_values.view();
    proc_cont.contacts = &qa_contacts;
    proc_cont.box = &box;
    proc_cont.r = r.view();
    proc_cont.free_pairs = &qa_free_pairs;

    auto &update = ker.update_qa_pairs;
    update.r = r.view();
    update.box = &box;
    update.nl = &nl;
    update.pairs = &qa_free_pairs;

    real cutoff = 0.0;
    for (int idx = 0; idx < qa::contact_type::NUM_TYPES; ++idx)
      cutoff = max(cutoff, sift.req_min_dist[idx]);
    nl_required = true;
    max_cutoff = max(max_cutoff, cutoff);

    if (params.out.enabled) {
      auto &report_qa = ker.report_qa_stuff;
      report_qa.sync_values = sync_values.view();
      report_qa.contacts = &qa_contacts;
      report_qa.process_cont = &ker.process_qa_contacts;
      hooks.push_back(&report_qa);
    }
  }
}

void simulation::setup_pid() {
  if (params.pid.enabled) {
    auto &eval = ker.eval_pid_forces;

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

    ss_ljs =
        nitro::vector<sink_lj>(amino_acid::NUM_TYPES * amino_acid::NUM_TYPES);
    for (auto const &aa1 : amino_acid::all()) {
      auto idx1 = (uint16_t)(uint8_t)aa1;
      for (auto const &aa2 : amino_acid::all()) {
        auto idx2 = (uint16_t)(uint8_t)aa2;
        auto ss_idx = idx1 * (uint16_t)amino_acid::NUM_TYPES + idx2;
        ss_ljs[ss_idx] = params.lj_variants.ss[{aa1, aa2}];
      }
    }
    eval.ss_ljs = ss_ljs.view();

    eval.r = r.view();
    eval.box = &box;
    eval.bundles = &pid_bundles;
    eval.prev = prev.view();
    eval.next = next.view();

    auto &update = ker.update_pid_bundles;
    update.r = r.view();
    update.prev = prev.view();
    update.next = next.view();
    update.atype = atype.view();
    update.box = &box;
    update.nl = &nl;
    update.bundles = &pid_bundles;

    real cutoff = 0.0;
    cutoff = max(cutoff, params.lj_variants.bb.cutoff());
    cutoff = max(cutoff, params.lj_variants.bs.cutoff());
    cutoff = max(cutoff, params.lj_variants.sb.cutoff());
    for (auto const &aa1 : amino_acid::all()) {
      for (auto const &aa2 : amino_acid::all()) {
        cutoff = max(cutoff, params.lj_variants.ss[{aa1, aa2}].cutoff());
      }
    }

    update.cutoff = cutoff;
    nl_required = true;
    max_cutoff = max(max_cutoff, cutoff);
  }
}

void simulation::setup_vafm() {
  if (params.vafm.enabled) {
    auto &eval = ker.eval_vafm_forces;
    eval.r = r.view();
    eval.t = &t;
    eval.afm_force.H1 = params.vafm.H1;
    eval.afm_force.H2 = params.vafm.H2;

    for (auto const &tip : params.vafm.afm_tips) {
      auto orig = tip.origin.value_or(r[tip.pulled_idx]);
      vafm_tips.emplace_back(tip.pulled_idx, orig, tip.v);
    }
    eval.afm_tips = vafm_tips.view();
  }
}

void simulation::setup_fafm() {
  if (params.fafm.enabled) {
    auto &eval = ker.eval_fafm_forces;

    for (auto const &tip : params.fafm.afm_tips) {
      fafm_tips.emplace_back(tip.pulled_idx, tip.F);
    }
    eval.afm_tips = fafm_tips.view();
  }
}

void simulation::setup_heur_ang() {
  if (params.heur_ang.enabled) {
    auto &eval = ker.eval_heur_ang_forces;
    eval.r = r.view();

    for (auto const &heur_pair : aa_heur_pair::all()) {
      for (int d = 0; d <= heur_ang::eval_forces::POLY_DEG; ++d) {
        eval.poly_coeffs[d][(uint8_t)heur_pair] =
            params.heur_ang.coeffs[heur_pair].poly[d];
      }
    }

    for (auto const &angle : model.angles) {
      if (!angle.theta.has_value()) {
        auto i1 = res_map[angle.res1], i2 = res_map[angle.res2],
             i3 = res_map[angle.res3];
        auto type = aa_heur_pair(atype[i2], atype[i3]);

        heur_angles.emplace_back(i1, i2, i3, type);
      }
    }
    eval.angles = heur_angles.view();
  }
}

void simulation::setup_heur_dih() {
  if (params.heur_dih.enabled) {
    auto &eval = ker.eval_heur_dih_forces;
    eval.r = r.view();

    for (auto const &heur_pair : aa_heur_pair::all()) {
      auto idx = (uint8_t)heur_pair;
      eval.coeffs.const_[idx] = params.heur_dih.coeffs[heur_pair].const_;
      eval.coeffs.sin[idx] = params.heur_dih.coeffs[heur_pair].sin;
      eval.coeffs.cos[idx] = params.heur_dih.coeffs[heur_pair].cos;
      eval.coeffs.sin2[idx] = params.heur_dih.coeffs[heur_pair].sin2;
      eval.coeffs.cos2[idx] = params.heur_dih.coeffs[heur_pair].cos2;
      eval.coeffs.sin_cos[idx] = params.heur_dih.coeffs[heur_pair].sin_cos;
    }

    for (auto const &dihedral : model.dihedrals) {
      if (!dihedral.phi.has_value()) {
        auto i1 = res_map[dihedral.res1], i2 = res_map[dihedral.res2],
             i3 = res_map[dihedral.res3], i4 = res_map[dihedral.res4];
        auto type = aa_heur_pair(atype[i2], atype[i3]);

        heur_dihedrals.emplace_back(i1, i2, i3, i4, type);
      }
    }
  }
}

void simulation::setup(int argc, char **argv) {
  parse_args(argc, argv);
  load_parameters();
  general_setup();
  load_model();
  compile_model();
  setup_dyn();
  setup_output();
  setup_langevin();
  setup_pbar();
  setup_nl();
  setup_chir();
  setup_tether();
  setup_nat_ang();
  setup_heur_ang();
  setup_nat_dih();
  setup_heur_dih();
  setup_pauli();
  setup_nat_cont();
  setup_dh();
  setup_qa();
  setup_pid();
  setup_vafm();
  setup_fafm();
}

simulation::thread::thread(simulation *simul) {
  this->simul = simul;
  params = simul->params;
  dyn = simul->dyn;
  gen = simul->gen;
  ker = simul->ker;

  ker.eval_chir_forces.F = dyn.F.view();
  ker.eval_chir_forces.V = &dyn.V;
  ker.eval_const_dh_forces.F = dyn.F.view();
  ker.eval_const_dh_forces.V = &dyn.V;
  ker.eval_rel_dh_forces.F = dyn.F.view();
  ker.eval_rel_dh_forces.V = &dyn.V;
  ker.eval_fafm_forces.F = dyn.F.view();
  ker.eval_heur_ang_forces.F = dyn.F.view();
  ker.eval_heur_ang_forces.V = &dyn.V;
  ker.eval_heur_dih_forces.F = dyn.F.view();
  ker.eval_heur_dih_forces.V = &dyn.V;
  ker.eval_nat_ang_forces.F = dyn.F.view();
  ker.eval_nat_ang_forces.V = &dyn.V;
  ker.eval_nat_cont_forces.F = dyn.F.view();
  ker.eval_nat_cont_forces.V = &dyn.V;
  ker.eval_cnd_forces.F = dyn.F.view();
  ker.eval_cnd_forces.V = &dyn.V;
  ker.eval_snd_forces.F = dyn.F.view();
  ker.eval_snd_forces.V = &dyn.V;
  ker.eval_pauli_forces.F = dyn.F.view();
  ker.eval_pauli_forces.V = &dyn.V;
  ker.eval_pid_forces.F = dyn.F.view();
  ker.eval_pid_forces.V = &dyn.V;
  ker.process_qa_contacts.F = dyn.F.view();
  ker.process_qa_contacts.V = &dyn.V;
  ker.eval_tether_forces.F = dyn.F.view();
  ker.eval_tether_forces.V = &dyn.V;
  ker.eval_vafm_forces.F = dyn.F.view();
  ker.lang_step.gen = &gen;
}

void simulation::thread::main_async() {
#pragma omp master
  pre_loop_init();

#pragma omp barrier
  real total_time = params.gen.total_time;
  while (simul->t < total_time) {
    async_part();
#pragma omp barrier

#pragma omp master
    sync_part();

#pragma omp barrier
  }
}

void simulation::thread::async_part() {
  dyn.reset();

  if (params.chir.enabled)
    ker.eval_chir_forces.omp_async();
  if (params.tether.enabled)
    ker.eval_tether_forces.omp_async();
  if (params.nat_ang.enabled)
    ker.eval_nat_ang_forces.omp_async();
  if (params.cnd.enabled)
    ker.eval_cnd_forces.omp_async();
  if (params.snd.enabled)
    ker.eval_snd_forces.omp_async();
  if (params.pauli.enabled)
    ker.eval_pauli_forces.omp_async();
  if (params.nat_cont.enabled)
    ker.eval_nat_cont_forces.omp_async();
  if (params.const_dh.enabled)
    ker.eval_const_dh_forces.omp_async();
  if (params.rel_dh.enabled)
    ker.eval_rel_dh_forces.omp_async();
  if (params.heur_ang.enabled)
    ker.eval_heur_ang_forces.omp_async();
  if (params.heur_dih.enabled)
    ker.eval_heur_dih_forces.omp_async();
  if (params.qa.enabled) {
    ker.prepare_nh.omp_async();
    ker.sift_qa_candidates.omp_async();
    ker.process_qa_contacts.omp_async();
  }
  if (params.pid.enabled)
    ker.eval_pid_forces.omp_async();
  if (params.vafm.enabled)
    ker.eval_vafm_forces.omp_async();
  if (params.fafm.enabled)
    ker.eval_fafm_forces.omp_async();

  dyn.omp_reduce(simul->dyn);
}

void simulation::thread::sync_part() {
  if (params.pbar.enabled)
    ker.render_pbar();

  if (params.out.enabled)
    ker.make_report();

  if (params.lang.enabled)
    ker.lang_step();
  simul->dyn.reset();
  if (params.qa.enabled)
    ker.sift_qa_candidates.omp_prep();

  if (simul->nl_required) {
    ker.nl_verify();
    if (simul->nl_invalid)
      on_nl_invalidation();
  }
}

void simulation::thread::pre_loop_init() {
#pragma omp master
  {
    simul->dyn.reset();
    if (params.qa.enabled)
      ker.sift_qa_candidates.omp_prep();

    if (simul->nl_required) {
      ker.nl_verify();
      if (simul->nl_invalid)
        on_nl_invalidation();
    }
  };
}

void simulation::thread::on_nl_invalidation() {
  if (params.nl.algorithm == nl::parameters::CELL) {
    ker.nl_cell();
  } else {
    ker.nl_legacy();
  }

  if (params.pauli.enabled)
    ker.update_pauli_pairs();
  if (params.nat_cont.enabled)
    ker.update_nat_contacts();
  if (params.const_dh.enabled || params.rel_dh.enabled)
    ker.update_dh_pairs();
  if (params.qa.enabled)
    ker.update_qa_pairs();
  if (params.pid.enabled)
    ker.update_pid_bundles();
}

int simulation::operator()(int argc, char **argv) {
  setup(argc, argv);
#pragma omp parallel default(none)
  {
    auto thr = thread(this);
    thr.main_async();
  }
  return EXIT_SUCCESS;
}