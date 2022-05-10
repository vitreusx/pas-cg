#include <Eigen/Dense>
#include <Eigen/SVD>
#include <cg/input/records.h>
#include <cg/output/make_report.h>
#include <fstream>
#include <iostream>

namespace cg::out {
static auto with_ext(std::filesystem::path const &prefix,
                     std::string const &ext) {
  return prefix.parent_path() / (prefix.filename().string() + ext);
}

static std::ofstream open_file(std::filesystem::path const &path) {
  if (path.has_parent_path()) {
    auto par = path.parent_path();
    if (!exists(par))
      std::filesystem::create_directories(par);
  }

  if (is_directory(path))
    throw std::runtime_error("file path \"" + path.string() +
                             "\" cannot be a directory");

  return std::ofstream(path);
}

void make_report::operator()() const {
  if (st->t >= rep->stats_last + stats_every) {
    add_cur_scalars();
    emit_out();
    add_map_data();
    emit_map();
    rep->stats_last = st->t;
  }

  if (st->t >= rep->struct_last + struct_every) {
    add_cur_snapshot();
    emit_pdb();
    rep->struct_last = st->t;
  }
}

real make_report::rmsd(nitro::const_view<vec3r> const &orig_r,
                       nitro::const_view<vec3r> const &cur_r,
                       const nitro::const_view<int> &indices) {
  auto n = indices.size();
  Eigen::MatrixX3d orig_R(n, 3), cur_R(n, 3);
  for (int idx = 0; idx < n; ++idx) {
    orig_R.row(idx) = convert<double>(orig_r[indices[idx]]).transpose();
    cur_R.row(idx) = convert<double>(cur_r[indices[idx]]).transpose();
  }

  Eigen::RowVector3d orig_center = orig_R.colwise().mean();
  orig_R.rowwise() -= orig_center;

  Eigen::RowVector3d cur_center = cur_R.colwise().mean();
  cur_R.rowwise() -= cur_center;

  Eigen::Matrix3d H = orig_R.transpose() * cur_R;
  auto svd = Eigen::JacobiSVD<decltype(H)>(H, Eigen::ComputeFullU |
                                                  Eigen::ComputeFullV);
  auto svdU = svd.matrixU(), svdV = svd.matrixV();
  auto d = sign((svdV.transpose() * svdU).determinant());
  auto R = svdV * Eigen::Vector3d(1, 1, d).asDiagonal() * svdU.transpose();

  auto rmsd = sqrt((cur_R * R - orig_R).rowwise().squaredNorm().sum());
  return (real)rmsd;
}

real make_report::kinetic_energy() const {
  real K = 0.0;
  for (int idx = 0; idx < st->r.size(); ++idx)
    K += (real)0.5 * st->comp_aa_data.mass[(uint8_t)st->atype[idx]] *
         norm_squared(st->v[idx]);
  return K;
}

real make_report::gyration_radius(const nitro::const_view<vec3r> &r,
                                  const nitro::const_view<int> &indices) {
  auto mean_r = vec3r::Zero();
  for (auto const &idx : indices)
    mean_r += r[idx];
  mean_r /= indices.size();

  real dr = 0.0;
  for (auto const &idx : indices)
    dr += norm_squared(r[idx] - mean_r);
  return sqrt(dr / (float)indices.size());
}

real make_report::asphericity(const nitro::const_view<vec3r> &r,
                              const nitro::const_view<int> &indices) {
  Eigen::MatrixX3d positions(indices.size(), 3);
  for (int idx = 0; idx < indices.size(); ++idx)
    positions.row(idx) = convert<double>(r[indices[idx]]).transpose();

  Eigen::RowVector3d center = positions.colwise().mean();
  positions.rowwise() -= center;

  Eigen::Vector3d lambda_ =
      Eigen::JacobiSVD<decltype(positions)>(positions).singularValues();
  lambda_ = lambda_.cwiseSqrt();
  auto lambda = vec3r(lambda_);

  return (real)1.5 * pow(lambda.x(), 2.0) - (real)0.5 * norm_squared(lambda);
}

void make_report::add_cur_scalars() const {
  auto &traj_div = rep->out_root.find<ioxx::sl4::div>(st->traj_idx);
  auto &scalars = traj_div.find<ioxx::sl4::table>("scalars");
  auto &row = scalars->append_row();

  nitro::vector<int> all_indices(st->r.size());
  for (int idx = 0; idx < st->r.size(); ++idx)
    all_indices[idx] = idx;

  auto full_Rg = gyration_radius(st->r, all_indices),
       full_W = asphericity(st->r, all_indices);
  auto full_rmsd = rmsd(st->orig_r, st->r, all_indices);
  auto K = kinetic_energy(), E = st->dyn.V + K;
  auto L = norm(st->r[st->r.size() - 1] - st->r[0]);

  row["TIME[tau]"] = quantity(st->t).in("tau");
  row["EPOT[eps]"] = quantity(st->dyn.V).in("eps");
  row["ETOT[eps]"] = quantity(E).in("eps");
  row["RG[A]"] = quantity(full_Rg).in("A");
  row["L[A]"] = quantity(L).in("A");
  row["RMSD[A]"] = quantity(full_rmsd).in("A");
  row["W[A**2]"] = quantity(full_W).in("A**2");

  if (nc) {
    int num_active_contacts = 0, num_active_ss = 0;
    for (auto const &cont : *nc->contacts) {
      if (nc->is_active(cont)) {
        ++num_active_contacts;
        if (cont.type() == nat_cont::type::SSBOND)
          ++num_active_ss;
      }
    }

    row["ICN"] = num_active_contacts;
    row["ICNss"] = num_active_ss;
  }
}

void make_report::emit_out() const {
  auto out_path = with_ext(prefix, ".out");
  auto out_of = open_file(out_path);
  out_of << rep->out_root;
}

void make_report::add_cur_snapshot() const {
  auto &snap = rep->snapshots.emplace_back();
  snap.t = st->t;
  snap.r = st->r;
  snap.model_box = st->box;
}

void make_report::emit_pdb() const {
  auto pdb_path = with_ext(prefix, ".pdb");
  auto pdb_of = open_file(pdb_path);

  std::stringstream ss{};

  bool first_line = true;
  for (int snap_idx = 0; snap_idx < (int)rep->snapshots.size(); ++snap_idx) {
    auto const &snap = rep->snapshots[snap_idx];

    auto chain_scalars_tab = ioxx::table::table();
    auto chain_indices = nitro::vector<int>();

    auto xmd_model = st->model;
    for (int idx = 0; idx < (int)xmd_model.residues.size(); ++idx) {
      auto *ref_res = st->model.residues[idx].get();
      auto res_idx = st->res_map.at(ref_res);
      xmd_model.residues[idx]->pos = snap.r[res_idx];
    }

    auto pdb_model = pdb_file(xmd_model).primary_model();
    pdb_model.model_serial = snap_idx + 1;

    for (auto const &[chain_id, chain] : pdb_model.chains) {
      auto chain_idx = chain_id - 'A';
      chain_indices.clear();
      for (int idx = 0; idx < st->num_res; ++idx)
        if (st->chain_idx[idx] == chain_idx)
          chain_indices.push_back(idx);

      auto chain_Rg = gyration_radius(snap.r, chain_indices);
      auto r_first = snap.r[st->chain_first[chain_idx]];
      auto r_last = snap.r[st->chain_last[chain_idx]];
      auto chain_L = norm(snap.model_box.wrap(r_first, r_last));
      auto chain_W = asphericity(snap.r, chain_indices);
      auto chain_N = chain_indices.size();

      auto &row = chain_scalars_tab.append_row();
      row["T[tau]"] = quantity(snap.t).in("tau");
      row["Rg[A]"] = quantity(chain_Rg).in("A");
      row["R_end_to_end[A]"] = quantity(chain_L).in("A");
      row["N"] = chain_N;
      row["W[A**2]"] = quantity(chain_W).in("A**2");
    }

    auto model_r = records::model();
    model_r.serial = pdb_model.model_serial;
    if (!first_line)
      pdb_of << '\n';
    pdb_of << model_r.write() << '\n';
    first_line = false;

    auto p = ioxx::table::sl4_parser();
    p.fit(chain_scalars_tab);

    ss = {};
    p.write(ss, chain_scalars_tab.cols, false);
    auto cols_str = ss.str();

    int row_idx = 0;
    for (auto const &[chain_id, chain] : pdb_model.chains) {
      for (auto const &rem : records::remark::create(1, cols_str))
        pdb_of << rem.write() << '\n';

      auto const &row = chain_scalars_tab.rows[row_idx];
      ss = {};
      p.write(ss, row, false);
      auto row_str = ss.str();

      for (auto const &rem : records::remark::create(1, row_str))
        pdb_of << rem.write() << '\n';

      pdb_of << chain;

      ++row_idx;
    }
  }
}

static std::string format_nc_type(nat_cont::type const &t) {
  switch (t) {
  case nat_cont::type::SSBOND:
    return "ssbond";
  case nat_cont::type::BACK_BACK:
    return "bb";
  case nat_cont::type::BACK_SIDE:
    return "bs";
  case nat_cont::type::SIDE_BACK:
    return "sb";
  case nat_cont::type::SIDE_SIDE:
    return "ss";
  case nat_cont::type::ANY:
    return "any";
  default:
    return "";
  }
}

static std::string format_qa_type(qa::contact_type const &t) {
  if (t == qa::contact_type::BACK_BACK())
    return "bb";
  else if (t == qa::contact_type::BACK_SIDE())
    return "bs";
  else if (t == qa::contact_type::SIDE_BACK())
    return "sb";
  else
    return "ss";
}

static std::string format_qa_status(qa::contact_status const &s) {
  switch (s) {
  case qa::contact_status::FORMING_OR_FORMED:
    return "formed/forming";
  case qa::contact_status::BREAKING:
    return "breaking";
  default:
    return "";
  }
}

static auto angle_value(cg::vec3r r1, cg::vec3r r2, cg::vec3r r3) {
  auto r12_u = unit(r2 - r1), r23_u = unit(r3 - r2);
  return acos(dot(r12_u, r23_u));
}

static auto dihedral_value(cg::vec3r r1, cg::vec3r r2, cg::vec3r r3,
                           cg::vec3r r4) {
  auto r12 = r2 - r1, r23 = r3 - r2, r34 = r4 - r3;
  auto x12_23 = cross(r12, r23), x23_34 = cross(r23, r34);
  auto x12_23_u = unit(x12_23), x23_34_u = unit(x23_34);
  auto phi = acos(dot(x12_23_u, x23_34_u));
  if (dot(x12_23, r34) < 0.0f)
    phi = -phi;
  return phi;
}

void make_report::add_map_data() const {
  auto &traj_div = rep->map_root.find<ioxx::sl4::div>(st->traj_idx);
  auto &cur_div = traj_div.add<ioxx::sl4::div>();

  cur_div.add<ioxx::sl4::comment>("T = ", st->t);

  if (nc) {
    cur_div.add<ioxx::sl4::comment>("native contacts");
    auto &num_comment = cur_div.add<ioxx::sl4::comment>("n = ", 0);
    auto &tab = cur_div.add<ioxx::sl4::table>();

    int num_active_contacts = 0;
    for (auto const &cont : *nc->contacts) {
      if (nc->is_active(cont)) {
        ++num_active_contacts;

        auto &nc_row = tab->append_row();
        nc_row["i1"] = cont.i1();
        nc_row["i2"] = cont.i2();
        nc_row["type"] = format_nc_type(cont.type());

        auto dist = norm(st->box.wrap(st->r[cont.i1()], st->r[cont.i2()]));
        nc_row["dist[A]"] = quantity(dist).in("A");
      }
    }

    num_comment = ioxx::sl4::comment("n = ", num_active_contacts);
  }

  if (qa) {
    cur_div.add<ioxx::sl4::comment>("qa");
    auto &num_comment = cur_div.add<ioxx::sl4::comment>("n = ", 0);
    auto &tab = cur_div.add<ioxx::sl4::table>();

    int num_active_contacts = 0;
    for (auto const &entry : *qa->contacts) {
      if (entry.is_vacant())
        continue;

      auto const &cont = entry.item();
      if (cont.status() == qa::FORMING_OR_FORMED) {
        ++num_active_contacts;

        auto &qa_row = tab->append_row();
        qa_row["i1"] = cont.i1();
        qa_row["i2"] = cont.i2();
        qa_row["type"] = format_qa_type(cont.type());
        qa_row["status"] = format_qa_status(cont.status());

        auto dist = norm(st->box.wrap(st->r[cont.i1()], st->r[cont.i2()]));
        qa_row["dist[A]"] = quantity(dist).in("A");
      }
    }

    num_comment = ioxx::sl4::comment("n = ", num_active_contacts);
  }

  if (pid) {
    cur_div.add<ioxx::sl4::comment>("pid");
    auto &num_comment = cur_div.add<ioxx::sl4::comment>("n = ", 0);
    auto &tab = cur_div.add<ioxx::sl4::table>();

    int num_active_contacts = 0;
    for (auto const &bundle : *pid->bundles) {
      if (pid->is_active(bundle)) {
        ++num_active_contacts;

        auto &pid_row = tab->append_row();
        pid_row["i1"] = bundle.i1();
        pid_row["i2"] = bundle.i2();

        auto dist = norm(st->box.wrap(st->r[bundle.i1()], st->r[bundle.i2()]));
        pid_row["dist[A]"] = quantity(dist).in("A");
      }
    }

    num_comment = ioxx::sl4::comment("n = ", num_active_contacts);
  }

  cur_div.add<ioxx::sl4::comment>("angles");
  auto &angles_tab = cur_div.add<ioxx::sl4::table>();
  for (auto const &angle : st->model.angles) {
    auto &row = angles_tab->append_row();

    auto i1 = st->res_map.at(angle.res1), i2 = st->res_map.at(angle.res2),
         i3 = st->res_map.at(angle.res3);
    row["i1"] = i1;
    row["i2"] = i2;
    row["i3"] = i3;

    auto r1 = st->r[i1], r2 = st->r[i2], r3 = st->r[i3];
    row["theta"] = angle_value(r1, r2, r3);

    if (angle.theta.has_value())
      row["native"] = angle.theta.value();
  }

  cur_div.add<ioxx::sl4::comment>("dihedrals");
  auto &dih_tab = cur_div.add<ioxx::sl4::table>();
  for (auto const &dih : st->model.dihedrals) {
    auto &row = dih_tab->append_row();

    auto i1 = st->res_map.at(dih.res1), i2 = st->res_map.at(dih.res2),
         i3 = st->res_map.at(dih.res3), i4 = st->res_map.at(dih.res4);
    row["i1"] = i1;
    row["i2"] = i2;
    row["i3"] = i3;
    row["i4"] = i4;

    auto r1 = st->r[i1], r2 = st->r[i2], r3 = st->r[i3], r4 = st->r[i4];
    row["phi"] = dihedral_value(r1, r2, r3, r4);

    if (dih.phi.has_value())
      row["native"] = dih.phi.value();
  }
}

void make_report::emit_map() const {
  auto map_path = with_ext(prefix, ".map");
  auto map_of = open_file(map_path);
  map_of << rep->map_root;
}
} // namespace cg::out
