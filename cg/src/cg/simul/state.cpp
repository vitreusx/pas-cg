#include <cg/simul/state.h>
namespace cg::simul {

state::state(ioxx::xyaml::node const &raw_params) {
  this->raw_params = raw_params.flatten();
  raw_params >> params;
}

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

static void allclose(vect::vector<vec3r> const &u,
                     vect::vector<vec3r> const &v) {
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
  traj_idx = 0;
  gen = rand_gen(params.gen.seed);
  load_model();
  simul_setup_output();
}

void state::simul_setup_output() {
  rep = out::report();
  out_enabled = params.out.enabled;
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

    if (source.ignore_cryst1)
      file.cryst1 = vec3<double>::Zero();

    orig_model = file.to_model();

    bond = 0.0;
    for (auto const &teth : orig_model.tethers)
      bond += teth.length.value();
    bond /= orig_model.tethers.size();

    if (params.input.morph_into_saw.has_value())
      params.input.morph_into_saw->bond_distance = bond;
  } else {
    auto &file = std::get<seq_file>(model_file_v);
    orig_model = std::move(file.model);
  }
  num_res = (int)orig_model.residues.size();

  if (!params.input.load_structure)
    orig_model.remove_native_structure();
}

void state::traj_init() {
  step_idx = 0;

  morph_model();
  compile_model();
  setup_dyn();
  setup_output();
  setup_langevin();
  setup_pbar();
  setup_nl();
  setup_local_rep();
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

  wall_type[X] = params.sbox.walls.x_axis;
  wall_type[Y] = params.sbox.walls.y_axis;
  wall_type[Z] = params.sbox.walls.z_axis;
  setup_walls();

  lj_walls_enabled = harmonic_walls_enabled = afm_enabled = false;
}

void state::morph_model() {
  model = orig_model;
  std::sort(
      model.contacts.begin(), model.contacts.end(),
      [](auto const &x, auto const &y) -> auto{
        auto p1 = std::make_pair(x.res1->seq_idx, x.res2->seq_idx);
        if (p1.first >= p1.second)
          std::swap(p1.first, p1.second);

        auto p2 = std::make_pair(y.res1->seq_idx, y.res2->seq_idx);
        if (p2.first >= p2.second)
          std::swap(p2.first, p2.second);

        return p1 < p2;
      });

  if (params.input.morph_into_saw.has_value()) {
    auto &saw_p = params.input.morph_into_saw.value();
    if (saw_p.perform) {
      model.morph_into_saw(gen, saw_p);
      if (std::holds_alternative<seq_file>(params.input.source))
        orig_model = model;
    }
  }

  //  if (params.input.morph_into_line.has_value()) {
  //    model.morph_into_line(params.input.morph_into_line.value());
  //  }
}

void state::compile_model() {
  int res_idx = 0;
  res_map = {};
  for (auto const &res : model.residues)
    res_map[res.get()] = res_idx++;

  r = vect::vector<vec3r>(num_res);
  for (auto const &res : model.residues)
    r[res_map.at(res.get())] = res->pos;

  decltype(res_map) orig_res_map;
  res_idx = 0;
  for (auto const &res : orig_model.residues)
    orig_res_map[res.get()] = res_idx++;

  orig_r = vect::vector<vec3r>(num_res);
  for (auto const &res : orig_model.residues)
    orig_r[orig_res_map.at(res.get())] = res->pos;

  atype = vect::vector<amino_acid>(num_res);
  for (auto const &res : model.residues)
    atype[res_map.at(res.get())] = res->type;

  comp_aa_data = compiled_aa_data(params.aa_data);
  if (params.input.normalize_mass) {
    for (auto &x : comp_aa_data.mass)
      x = 1.0;
  }

  auto init_p = params.sbox.init_size;
  if (init_p.type == "from CRYST1 record") {
    if (!model.cryst1.has_value())
      throw std::runtime_error("CRYST1 record has not been found");

    box.min = vec3r::Zero();
    box.max = model.cryst1.value();
  } else if (init_p.type == "sufficient") {
    box = sbox::box<real>();
    auto suf_p = init_p.sufficient;
    real pad = suf_p.pad * bond;
    real side = 0;

    if (suf_p.max_density.has_value()) {
      real vol = num_res / suf_p.max_density.value();
      side = cbrt(vol) / (real)2.0;
      box.max = {side, side, side};
      box.min = -box.max;
    };

    if (suf_p.cubic || params.sbox.squeezing.perform) {
      for (int idx = 0; idx < num_res; ++idx) {
        auto p = r[idx];
        side = max(side, abs(p.x()) + pad);
        side = max(side, abs(p.y()) + pad);
        side = max(side, abs(p.z()) + pad);
      }
      box.max = {side, side, side};
      box.min = -box.max;
    } else {
      for (int idx = 1; idx < num_res; ++idx) {
        auto p = r[idx];
        box.min.x() = min(box.min.x(), p.x() - pad);
        box.min.y() = min(box.min.y(), p.y() - pad);
        box.min.z() = min(box.min.z(), p.z() - pad);
        box.max.x() = max(box.max.x(), p.x() + pad);
        box.max.y() = max(box.max.y(), p.y() + pad);
        box.max.z() = max(box.max.z(), p.z() + pad);
      }
    }
  } else if (init_p.type == "infinite") {
    auto void_x = params.sbox.walls.x_axis == "void";
    auto void_y = params.sbox.walls.y_axis == "void";
    auto void_z = params.sbox.walls.z_axis == "void";
    if (!void_x || !void_y || !void_z)
      throw std::runtime_error("for an infinite box, walls must be void");
  } else if (init_p.type == "keep from SAW") {
    if (!model.saw_box.has_value())
      throw std::runtime_error("getting initial box size from SAW requires SAW "
                               "with periodic boundary conditions to be run");
    box.min = model.saw_box->min;
    box.max = model.saw_box->max;
  }

  prev = next = chain_idx = seq_idx = vect::vector<int>(num_res, -1);
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

  chain_first = chain_last = vect::vector<int>((int)model.chains.size());
  for (auto const &chain : model.chains) {
    chain_first[chain->chain_idx] = res_map.at(chain->residues.front());
    chain_last[chain->chain_idx] = res_map.at(chain->residues.back());
  }
}

void state::setup_dyn() {
  t = 0;
  dyn = dynamics(num_res);
}

void state::setup_output() {
  out_enabled = params.out.enabled;
  if (out_enabled)
    rep.traj_init(traj_idx);

  ckpt_enabled = params.ckpt.enabled;
  if (ckpt_enabled)
    ckpt_last_t = std::numeric_limits<real>::lowest();

  dump_enabled = params.gen.dump_data;
  if (dump_enabled)
    std::ofstream raw_data_txt("raw_data.txt");
}

void state::setup_langevin() {
  lang_enabled = params.lang.enabled;

  if (lang_enabled) {
    auto &mass = comp_aa_data.mass;

    mass_inv = vect::vector<real>(mass.size());
    for (int aa_idx = 0; aa_idx < mass_inv.size(); ++aa_idx)
      mass_inv[aa_idx] = 1.0 / mass[aa_idx];

    mass_rsqrt = vect::vector<real>(mass.size());
    for (int aa_idx = 0; aa_idx < mass_rsqrt.size(); ++aa_idx)
      mass_rsqrt[aa_idx] = sqrt(mass_inv[aa_idx]);

    v = vect::vector<vec3r>(num_res, vec3r::Zero());

    y0 = y1 = y2 = y3 = y4 = y5 = vect::vector<vec3sr>(num_res);
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

    {
      vec3sr sum_vel;
      for (int idx = 0; idx < num_res; ++idx) {
        auto xx = 2.0 * (gen.uniform<real>() - 0.5);
        auto yy = 2.0 * (gen.uniform<real>() - 0.5);
        auto zz = 2.0 * (gen.uniform<real>() - 0.5);
        auto xyz = 1.0 / sqrt(xx * xx + yy * yy + zz * zz);

        y1[idx] = vec3sr(xx, yy, zz) * xyz;
        sum_vel += y1[idx];
      }

      auto x = 0.0;
      for (int idx = 0; idx < num_res; ++idx) {
        y1[idx] -= sum_vel / num_res;
        //      x += norm_squared(y1[idx]);
        decltype(auto) v = y1[idx];
        x = x + v.x() * v.x() + v.y() * v.y() + v.z() * v.z();
      }

      auto dt = params.lang.dt;
      auto delsq = dt * dt;
      auto aheat = delsq * 3.0 * num_res * temperature;
      auto heat = sqrt(aheat / x);
      for (int idx = 0; idx < num_res; ++idx) {
        y1[idx] *= heat;
      }
    }

    noise = vect::vector<vec3r>(num_res);
  }
}

void state::setup_pbar() {
  pbar_enabled = params.pbar.enabled && !dump_enabled;

  if (pbar_enabled) {
    pbar_first_time = true;
  }
}

void state::setup_nl() {
  nl_invalid = true;
  nl_required = false;
  verify_first_time = true;
  total_disp = 0.0;
  nl = nl::data();
  nl.orig_r = vect::vector<vec3r>(num_res);
}

void state::setup_local_rep() {
  lrep_enabled = params.lrep.enabled;

  if (lrep_enabled) {
    local_rep_pairs = {};
    for (auto const &triple : model.angles) {
      auto i1 = res_map[triple.res1], i3 = res_map[triple.res3];
      local_rep_pairs.emplace_back(i1, i3);
    }
  }
}

void state::setup_chir() {
  chir_enabled = params.chir.enabled;

  if (chir_enabled) {
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
  tether_enabled = params.tether.enabled;

  if (tether_enabled) {
    tether_pairs = {};
    for (auto const &tether : model.tethers) {
      auto i1 = res_map[tether.res1], i2 = res_map[tether.res2];
      //      auto nat_dist =
      //      (real)tether.length.value_or(params.tether.def_length);
      auto nat_dist = norm(pbc.wrap<vec3r>(orig_r[i1], orig_r[i2]));
      tether_pairs.emplace_back(i1, i2, nat_dist);
    }
  }
}

void state::setup_nat_ang() {
  nat_ang_enabled = params.angles.nat_ang.enabled;

  if (nat_ang_enabled) {
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
  heur_ang_enabled = params.angles.heur_ang.enabled;

  if (heur_ang_enabled) {
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
  nat_dih_enabled = params.angles.nat_dih.enabled;

  if (nat_dih_enabled) {
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
  heur_dih_enabled = params.angles.heur_dih.enabled;

  if (heur_dih_enabled) {
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
  pauli_enabled = params.pauli.enabled && !params.qa.enabled;

  if (pauli_enabled) {
    nl_required = true;
  }
}

void state::setup_nat_cont() {
  nat_cont_enabled = !model.contacts.empty() && params.nat_cont.enabled;

  if (nat_cont_enabled) {
    all_native_contacts = {};
    nat_cont_excl = {};
    num_changed = 0;
    for (int idx = 0; idx < (int)model.contacts.size(); ++idx) {
      auto const &cont = model.contacts[idx];
      auto i1 = res_map[cont.res1], i2 = res_map[cont.res2];
      if (i1 > i2)
        std::swap(i1, i2);

      auto nat_dist = (real)cont.length;
      bool active = norm(r[i2] - r[i1]) <
                    params.gen.counting_factor * C216_INV * nat_dist;
      real change_t = -1.0;
      all_native_contacts.emplace_back(i1, i2, nat_dist, cont.type, active,
                                       change_t, idx);
      nat_cont_excl.emplace_back(i1, i2);
    }
  }
}

void state::setup_dh() {
  dh_enabled = params.dh.enabled;

  if (dh_enabled) {
    nl_required = true;
  }
}

void state::setup_qa() {
  qa_enabled = params.qa.enabled;
  if (qa_enabled) {
    sync_values = vect::vector<sync_data>(num_res);
    for (int res_idx = 0; res_idx < num_res; ++res_idx) {
      auto lim = comp_aa_data.base_sync[(uint8_t)atype[res_idx]];
      sync_values[res_idx] = sync_data(
          lim.back(), lim.side_all(), lim.side_polar(), lim.side_hydrophobic());
    }

    n = h = vect::vector<vec3r>(num_res);
    num_qa_contacts = 0;

    if (params.qa.disulfide.has_value()) {
      neigh_count = vect::vector<int>(num_res, 0);
      part_of_ssbond = vect::vector<bool>(num_res, false);

      cys_indices = vect::vector<int>(num_res);
      for (int idx = 0; idx < num_res; ++idx) {
        if (atype[idx] == amino_acid(aa_code::CYS))
          cys_indices.push_back(idx);
      }
    }

    nl_required = true;
  }
}

void state::setup_pid() {
  pid_enabled = params.pid.enabled;

  if (pid_enabled) {
    nl_required = true;
  }
}

void state::setup_afm() {
  if (afm_enabled) {
    auto i0 = chain_first[0], i1 = chain_last[0];

    auto avg = moving_avg<vec3r, real>(params.afm.avg_stats_over);
    auto perp_avg = moving_avg<real, real>(params.afm.avg_stats_over);

    vel_afm_tips.emplace_back(i0, r[i0], t, vec3r::Zero(), avg, perp_avg);

    auto pull_dir = unit(r[i1] - r[i0]);
    auto type = params.afm.tip_params.type;
    if (type == "const velocity" || cur_phase == phase::PULL_RELEASE) {
      auto pull_vel = params.afm.tip_params.vel_afm.vel;
      vel_afm_tips.emplace_back(i1, r[i1], t, pull_vel * pull_dir, avg,
                                perp_avg);
    } else if (type == "const force") {
      auto pull_force = params.afm.tip_params.force_afm.force;
      force_afm_tips.emplace_back(i1, pull_force * pull_dir, avg);
    }
  }
}

void state::setup_walls() {
  harmonic_walls.clear();
  harmonic_conns.clear();
  dyn.harmonic_wall_F.clear();

  solid_walls.clear();
  dyn.solid_wall_F.clear();

  ljw_removed.clear();
  ljw_conns.clear();
  ljw_candidates.clear();
  lj_walls.clear();
  dyn.lj_wall_F.clear();

  is_connected_to_wall = vect::vector<bool>(num_res, false);
  harmonic_walls_enabled = solid_walls_enabled = lj_walls_enabled = false;
  walls = vect::vector<wall::gen_wall *>(6, nullptr);

  for (auto const &axis : {X, Y, Z}) {
    vec3r normal;
    side neg_side, pos_side;

    switch (axis) {
    case X:
      normal = vec3r::UnitX();
      neg_side = NEG_X;
      pos_side = POS_X;
      break;
    case Y:
      normal = vec3r::UnitY();
      neg_side = NEG_Y;
      pos_side = POS_Y;
      break;
    case Z:
      normal = vec3r::UnitZ();
      neg_side = NEG_Z;
      pos_side = POS_Z;
      break;
    }

    auto neg_axis_plane = cg::plane<real>(vec3r::Zero(), normal);
    auto pos_axis_plane = cg::plane<real>(vec3r::Zero(), -normal);

    int auto_limit =
        (num_res -
         nearbyint(
             pow(cbrt(num_res) - params.sbox.walls.threshold *
                                     cbrt(params.sbox.squeezing.target_density),
                 3.0))) /
        3;
    auto_limit /= 2;

    if (wall_type[axis] == "solid") {
      decltype(auto) neg_force = dyn.solid_wall_F.emplace_back();
      decltype(auto) neg_wall = solid_walls.emplace_back(
          neg_axis_plane, vect::iterator<vec3r>(neg_force),
          params.sbox.avg_forces_over);
      walls[neg_side] = &neg_wall;

      decltype(auto) pos_force = dyn.solid_wall_F.emplace_back();
      decltype(auto) pos_wall = solid_walls.emplace_back(
          pos_axis_plane, vect::iterator<vec3r>(pos_force),
          params.sbox.avg_forces_over);
      walls[pos_side] = &pos_wall;

      solid_walls_enabled = true;
    } else if (wall_type[axis] == "harmonic") {
      auto limit = params.sbox.walls.harmonic_wall.limit.value_or(auto_limit);

      decltype(auto) neg_force = dyn.harmonic_wall_F.emplace_back();
      decltype(auto) neg_wall = harmonic_walls.emplace_back(
          neg_axis_plane, vect::iterator<vec3r>(neg_force),
          params.sbox.avg_forces_over, limit);
      walls[neg_side] = &neg_wall;

      decltype(auto) pos_force = dyn.harmonic_wall_F.emplace_back();
      decltype(auto) pos_wall = harmonic_walls.emplace_back(
          pos_axis_plane, vect::iterator<vec3r>(pos_force),
          params.sbox.avg_forces_over, limit);
      walls[pos_side] = &pos_wall;

      harmonic_walls_enabled = true;

      vect::vector<std::pair<real, int>> dist(num_res);
      for (int wall_idx = 0; wall_idx < harmonic_walls.size(); ++wall_idx) {
        decltype(auto) wall = harmonic_walls[wall_idx];
        for (int res_idx = 0; res_idx < num_res; ++res_idx) {
          dist[res_idx] = std::make_pair(wall.plane.dist(r[res_idx]), res_idx);
        }
        std::sort(dist.begin(), dist.end());

        for (int idx = 0, conn = 0; conn < limit && idx < num_res; ++idx) {
          auto res_idx = dist[idx].second;
          if (!is_connected_to_wall[res_idx]) {
            auto proj = wall.plane.projection(r[res_idx]);
            auto offset = proj - wall.plane.origin();
            real saturation = 0.0;
            harmonic_conns.emplace_back(wall_idx, res_idx, offset, saturation);
            is_connected_to_wall[res_idx] = true;
            ++conn;
          }
        }
      }
    } else if (wall_type[axis] == "lj") {
      auto limit = params.sbox.walls.lj_wall.limit.value_or(auto_limit);

      decltype(auto) neg_force = dyn.lj_wall_F.emplace_back();
      decltype(auto) neg_wall = lj_walls.emplace_back(
          neg_axis_plane, vect::iterator<vec3r>(neg_force),
          params.sbox.avg_forces_over, limit);
      walls.push_back(&neg_wall);
      walls[neg_side] = &neg_wall;

      decltype(auto) pos_force = dyn.lj_wall_F.emplace_back();
      decltype(auto) pos_wall = lj_walls.emplace_back(
          pos_axis_plane, vect::iterator<vec3r>(pos_force),
          params.sbox.avg_forces_over, limit);
      walls[pos_side] = &pos_wall;

      lj_walls_enabled = true;
    }

    pbc_on[axis] = (wall_type[axis] == "periodic");
  }

  if (walls[NEG_Z] && walls[POS_Z]) {
    avg_z_force = moving_avg<real, real>(params.sbox.avg_forces_over);
  }

  reinit_wall_values();
}

void state::adjust_wall_pos(vec3r size_change, vec3r translation) {
  box.min += -size_change + translation;
  box.max += size_change + translation;
  reinit_wall_values();
}

void state::reinit_wall_values() {
  vec3r center = box.center(), ext = box.extent();

  for (auto side : {NEG_X, POS_X, NEG_Y, POS_Y, NEG_Z, POS_Z}) {
    if (!walls[side])
      continue;

    auto prev_val = walls[side]->plane.origin();
    auto is_pos = side % 2 == 1;
    auto next_val =
        center + cg::cast(is_pos ? box.max - center : box.min - center,
                          walls[side]->plane.normal());
    walls[side]->plane.origin() = next_val;
    walls[side]->shift = next_val - prev_val;
  }

  vec3r pbc_cell;
  pbc_cell.x() = pbc_on[X] ? ext.x() : 0;
  pbc_cell.y() = pbc_on[Y] ? ext.y() : 0;
  pbc_cell.z() = pbc_on[Z] ? ext.z() : 0;
  pbc.set_cell(pbc_cell);
}

bool state::trajectory_should_end() const {
  if (cur_phase == phase::TRAJ_INIT || cur_phase == phase::SIMUL_INIT ||
      cur_phase == phase::SIMUL_END)
    return false;
  if (t >= params.gen.total_time)
    return true;
  if (params.nat_cont.unfolding_study.early_stopping) {
    if (num_changed == all_native_contacts.size())
      return true;
  }
  return false;
}

} // namespace cg::simul