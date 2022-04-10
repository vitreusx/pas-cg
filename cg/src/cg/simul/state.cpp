#include <cg/simul/state.h>
namespace cg::simul {

state::state(parameters &params) : params(params) {}

static void allclose(real u, real v) {
  if (abs(u) < 1e-8 && abs(v) < 1e-8)
    return;

  if (abs(v - u) > 1e-8 * abs(u))
    throw;
}

static void allclose(vec3r u, vec3r v) {
  allclose(u.x(), v.x());
  allclose(u.y(), v.y());
  allclose(u.z(), v.z());
}

static void allclose(nitro::vector<vec3r> const &u,
                     nitro::vector<vec3r> const &v) {
  if (u.size() != v.size())
    throw;

  for (int idx = 0; idx < u.size(); ++idx) {
    allclose(u[idx], v[idx]);
  }
}

void state::verify_equal(const state &other) const {
  allclose(dyn.V, other.dyn.V);
  allclose(dyn.F, other.dyn.F);
  allclose(r, other.r);
}

void state::simul_setup() {
  is_running = true;
  total_time = params.gen.total_time;
  equil_time = params.gen.equil_time;
  traj_idx = 0;

  if (params.gen.disable_all) {
    params.chir.enabled = false;
    params.tether.enabled = false;
    params.nat_ang.enabled = false;
    params.heur_ang.enabled = false;
    params.cnd.enabled = false;
    params.snd.enabled = false;
    params.heur_dih.enabled = false;
    params.pauli.enabled = false;
    params.nat_cont.enabled = false;
    params.const_dh.enabled = false;
    params.rel_dh.enabled = false;
    params.qa.enabled = false;
    params.pid.enabled = false;
  }

  gen = rand_gen(params.gen.seed);

  report = out::report_data();
  report.out_dir = params.out.output_dir;
  report.last_stats_t = std::numeric_limits<real>::lowest();
  report.last_files_t = report.last_stats_t;
  report.traj_idx = &traj_idx;
  report.simul_first_time = true;

  load_model();
}

void state::load_model() {
  auto &model_file_v = params.input.source;
  if (std::holds_alternative<input::parameters::pdb_source>(model_file_v)) {
    auto const &source = std::get<input::parameters::pdb_source>(model_file_v);
    auto file = source.file;
    if (source.deriv) {
      auto all_atoms = source.deriv == pdb_file::contact_deriv::FROM_ATOMS;
      file.add_contacts(params.aa_data, all_atoms);
    }
    orig_model = file.to_model();
  } else {
    auto &file = std::get<seq_file>(model_file_v);
    orig_model = std::move(file.model);
  }
  num_res = (int)orig_model.residues.size();
}

void state::traj_setup() {
  morph_model();
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
  setup_afm();
}

void state::finish_trajectory() {
  did_traj_setup = false;
  did_post_equil_setup = false;
  ++traj_idx;
  gen = gen.spawn();
}

void state::morph_model() {
  model = orig_model;

  if (params.input.morph_into_saw.has_value()) {
    auto &saw_p = params.input.morph_into_saw.value();
    if (saw_p.perform)
      model.morph_into_saw(gen, saw_p);
  }
  if (params.input.morph_into_line.has_value()) {
    model.morph_into_line(params.input.morph_into_line.value());
  }
}

void state::compile_model() {
  int res_idx = 0;
  res_map = {};
  for (auto const &res : model.residues)
    res_map[res.get()] = res_idx++;

  r = nitro::vector<vec3r>(num_res);
  for (auto const &res : model.residues)
    r[res_map.at(res.get())] = res->pos;

  decltype(res_map) orig_res_map;
  res_idx = 0;
  for (auto const &res : orig_model.residues)
    orig_res_map[res.get()] = res_idx++;

  orig_r = nitro::vector<vec3r>(num_res);
  for (auto const &res : orig_model.residues)
    orig_r[orig_res_map.at(res.get())] = res->pos;

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

  chain_first = chain_last = nitro::vector<int>((int)model.chains.size());
  for (auto const &chain : model.chains) {
    chain_first[chain->chain_idx] = res_map.at(chain->residues.front());
    chain_last[chain->chain_idx] = res_map.at(chain->residues.back());
  }
}

void state::setup_dyn() {
  t = 0;
  dyn = dynamics(num_res);
}

void state::setup_output() { report.on_new_trajectory(); }

void state::setup_langevin() {
  auto &mass = comp_aa_data.mass;

  mass_inv = nitro::vector<real>(mass.size());
  for (int aa_idx = 0; aa_idx < mass_inv.size(); ++aa_idx)
    mass_inv[aa_idx] = 1.0 / mass[aa_idx];

  mass_rsqrt = nitro::vector<real>(mass.size());
  for (int aa_idx = 0; aa_idx < mass_rsqrt.size(); ++aa_idx)
    mass_rsqrt[aa_idx] = sqrt(mass_inv[aa_idx]);

  v = nitro::vector<vec3r>(num_res, vec3r::Zero());

  y0 = y1 = y2 = y3 = y4 = y5 = nitro::vector<vec3sr>(num_res);
  for (int idx = 0; idx < num_res; ++idx)
    y0[idx] = r[idx];

  true_t = t;

  if (std::holds_alternative<quantity>(params.lang.temperature)) {
    auto temp = std::get<quantity>(params.lang.temperature);
    temperature = temp;
  } else {
    auto [temp_start_q, temp_end_q] =
        std::get<lang::parameters::quantity_range>(params.lang.temperature);
    real temp_start = temp_start_q, temp_end = temp_end_q;

    auto w = (real)traj_idx / (real)(params.gen.num_of_traj - 1);
    temperature = temp_start * w + temp_end * ((real)1.0 - w);
  }
}

void state::setup_pbar() {
  if (params.pbar.enabled) {
    pbar_first_time = true;
  }
}

void state::setup_nl() {
  nl_invalid = true;
  nl_required = false;
  verify_first_time = true;
  total_disp = 0.0;
  max_cutoff = 0.0;
  nl = nl::data();
  nl.orig_r = nitro::vector<vec3r>(num_res);

  switch (params.nl.algorithm) {
  case nl::parameters::LEGACY:
    break;
  case nl::parameters::CELL:
    res_cell_idx = reordered_idx = nitro::vector<int>(num_res);
    break;
  }
}

void state::setup_chir() {
  if (params.chir.enabled) {
    chir_quads = {};
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
  }
}

void state::setup_tether() {
  if (params.tether.enabled) {
    tether_pairs = {};
    for (auto const &tether : model.tethers) {
      auto i1 = res_map[tether.res1], i2 = res_map[tether.res2];
      auto nat_dist = (real)tether.length.value_or(params.tether.def_length);
      tether_pairs.emplace_back(i1, i2, nat_dist);
    }
  }
}

void state::setup_nat_ang() {
  if (params.nat_ang.enabled) {
    native_angles = {};
    for (auto const &angle : model.angles) {
      if (angle.theta.has_value()) {
        auto i1 = res_map[angle.res1], i2 = res_map[angle.res2],
             i3 = res_map[angle.res3];
        auto nat_theta = (real)angle.theta.value();

        native_angles.emplace_back(i1, i2, i3, nat_theta);
      }
    }
  }
}

void state::setup_heur_ang() {
  if (params.heur_ang.enabled) {
    heur_angles = {};
    for (auto const &angle : model.angles) {
      if (!angle.theta.has_value()) {
        auto i1 = res_map[angle.res1], i2 = res_map[angle.res2],
             i3 = res_map[angle.res3];
        auto type = aa_heur_pair(atype[i2], atype[i3]);

        heur_angles.emplace_back(i1, i2, i3, type);
      }
    }
  }
}

void state::setup_nat_dih() {
  if (params.cnd.enabled || params.snd.enabled) {
    native_dihedrals = {};
    for (auto const &dihedral : model.dihedrals) {
      if (dihedral.phi.has_value()) {
        auto i1 = res_map[dihedral.res1], i2 = res_map[dihedral.res2],
             i3 = res_map[dihedral.res3], i4 = res_map[dihedral.res4];
        auto nat_phi = (real)dihedral.phi.value();

        native_dihedrals.emplace_back(i1, i2, i3, i4, nat_phi);
      }
    }
  }
}

void state::setup_heur_dih() {
  if (params.heur_dih.enabled) {
    heur_dihedrals = {};
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

void state::setup_pauli() {
  if (params.pauli.enabled) {
    nl_required = true;
    pauli_cutoff = (real)params.pauli.r_excl;
    max_cutoff = max(max_cutoff, pauli_cutoff);
  }
}

void state::setup_nat_cont() {
  if (model.contacts.empty())
    params.nat_cont.enabled = false;

  if (params.nat_cont.enabled) {
    all_native_contacts = {};
    nat_cont_excl = {};
    for (int idx = 0; idx < (int)model.contacts.size(); ++idx) {
      auto const &cont = model.contacts[idx];
      auto i1 = res_map[cont.res1], i2 = res_map[cont.res2];
      auto nat_dist = (real)cont.length;
      bool formed = false;
      real formation_t = 0.0;
      all_native_contacts.emplace_back(i1, i2, nat_dist, cont.type, formed,
                                       formation_t, idx);
      nat_cont_excl.emplace_back(i1, i2);
    }

    nat_cont_cutoff = 0.0;
    for (int cont_idx = 0; cont_idx < all_native_contacts.size(); ++cont_idx) {
      auto cont_dist = all_native_contacts[cont_idx].nat_dist();
      nat_cont_cutoff = max(nat_cont_cutoff, lj::compute_cutoff(cont_dist));
    }
    nl_required = true;
    max_cutoff = max(max_cutoff, nat_cont_cutoff);
  }
}

void state::setup_dh() {
  if (params.const_dh.enabled || params.rel_dh.enabled) {
    nl_required = true;

    real dh_cutoff = 0.0f;
    if (params.const_dh.enabled) {
      dh_cutoff = 2.0 * params.const_dh.screening_dist;
    } else {
      dh_cutoff = 2.0 * params.rel_dh.screening_dist;
    }

    max_cutoff = max(max_cutoff, dh_cutoff);
  }
}

void state::setup_qa() {
  if (params.qa.enabled) {
    sync_values = nitro::vector<sync_data>(num_res);
    for (int res_idx = 0; res_idx < num_res; ++res_idx) {
      auto lim = comp_aa_data.base_sync[(uint8_t)atype[res_idx]];
      sync_values[res_idx] = sync_data(
          lim.back(), lim.side_all(), lim.side_polar(), lim.side_hydrophobic());
    }

    n = h = nitro::vector<vec3r>(num_res);
    num_qa_contacts = 0;

    if (params.qa.disulfide.has_value()) {
      neigh_count = nitro::vector<int>(num_res, 0);
      part_of_ssbond = nitro::vector<bool>(num_res, false);

      cys_indices = nitro::vector<int>(num_res);
      for (int idx = 0; idx < num_res; ++idx) {
        if (atype[idx] == amino_acid(aa_code::CYS))
          cys_indices.push_back(idx);
      }
    }

    nl_required = true;
  }
}

void state::setup_pid() {
  if (params.pid.enabled) {
    ss_ljs =
        nitro::vector<sink_lj>(amino_acid::NUM_TYPES * amino_acid::NUM_TYPES);
    for (auto const &aa1 : amino_acid::all()) {
      auto idx1 = (uint16_t)(uint8_t)aa1;
      for (auto const &aa2 : amino_acid::all()) {
        auto idx2 = (uint16_t)(uint8_t)aa2;
        auto ss_idx = idx1 * (uint16_t)amino_acid::NUM_TYPES + idx2;
        ss_ljs[ss_idx] = params.lj_vars.ss[{aa1, aa2}];
      }
    }

    nl_required = true;
  }
}

void state::setup_afm() {
  if (!params.afm.tips.empty()) {
    afm_tips = afm::compiled_tips();

    for (auto const &tip : params.afm.tips) {
      if (std::holds_alternative<afm::parameters::single_res_t>(tip)) {
        auto const &tip_v = std::get<afm::parameters::single_res_t>(tip);

        if (tip_v.type == afm::parameters::tip_type::CONST_VEL) {
          auto orig = r[tip_v.res_idx];
          afm_tips.vel.emplace_back(tip_v.res_idx, orig, (vec3r)tip_v.dir);
        } else {
          afm_tips.force.emplace_back(tip_v.res_idx, (vec3r)tip_v.dir);
        }

      } else {
        auto const &tip_v = std::get<afm::parameters::pulled_apart_t>(tip);
        afm_tips.pulled_chains.push_back(tip_v.chain_idx);
        auto first = chain_first[tip_v.chain_idx],
             last = chain_last[tip_v.chain_idx];
        auto r_first = r[first], r_last = r[last];
        auto dir = unit(r_last - r_first);

        if (tip_v.type == afm::parameters::tip_type::CONST_VEL) {
          afm_tips.vel.emplace_back(first, r_first, -dir * tip_v.mag);
          afm_tips.vel.emplace_back(last, r_last, dir * tip_v.mag);
        } else {
          afm_tips.force.emplace_back(first, -dir * tip_v.mag);
          afm_tips.force.emplace_back(last, dir * tip_v.mag);
        }
      }
    }
  }
}

void state::post_equil_setup() {
  setup_afm();
  post_equil = true;
}

} // namespace cg::simul