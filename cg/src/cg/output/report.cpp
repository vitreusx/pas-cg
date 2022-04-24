#include <Eigen/Dense>
#include <Eigen/SVD>
#include <cg/output/report.h>
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

void report::simul_init() { model_serial = 1; }

void report::traj_init(int traj_idx) {
  stats_last = std::numeric_limits<real>::lowest();
  struct_last = std::numeric_limits<real>::lowest();

  auto &traj_out_div = out_root.add<ioxx::sl4::div>();
  traj_out_div.add<ioxx::sl4::comment>("trajectory ", traj_idx);
  traj_out_div.named_add<ioxx::sl4::table>("scalars");

  auto &traj_map_div = map_root.add<ioxx::sl4::div>();
  traj_map_div.add<ioxx::sl4::comment>("trajectory ", traj_idx);
}

void make_report::operator()() const {
  if (*t >= rep->stats_last + stats_every) {
    add_cur_scalars();
    emit_out();
    add_map_data();
    emit_map();
    rep->stats_last = *t;
  }

  if (*t >= rep->struct_last + struct_every) {
    add_cur_model();
    emit_pdb();
    rep->struct_last = *t;
  }
}

gyration::gyration(const nitro::const_view<vec3r> &r,
                   const nitro::const_view<int> &indices) {
  Eigen::MatrixX3d positions(indices.size(), 3);
  for (int idx = 0; idx < indices.size(); ++idx)
    positions.row(idx) = convert<double>(r[indices[idx]]).transpose();

  Eigen::RowVector3d center = positions.colwise().mean();
  positions.rowwise() -= center;

  Eigen::Vector3d lambda_ =
      Eigen::JacobiSVD<decltype(positions)>(positions).singularValues();
  lambda_ = lambda_.cwiseSqrt();
  lambda = vec3r(lambda_);
}

real gyration::radius() const { return norm(lambda); }

real gyration::asphericity() const {
  return (real)1.5 * pow(lambda.x(), 2.0) - (real)0.5 * norm_squared(lambda);
}

real make_report::rmsd(const nitro::const_view<int> &indices) const {
  auto n = indices.size();
  Eigen::MatrixX3d orig_R(n, 3), cur_R(n, 3);
  for (int idx = 0; idx < n; ++idx) {
    orig_R.row(idx) = convert<double>(orig_r[indices[idx]]).transpose();
    cur_R.row(idx) = convert<double>(r[indices[idx]]).transpose();
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
  for (int idx = 0; idx < r.size(); ++idx)
    K += (real)0.5 * mass[(uint8_t)atype[idx]] * norm_squared(v[idx]);
  return K;
}

real make_report::gyration_radius(const nitro::const_view<int> &indices) const {
  auto mean_r = vec3r::Zero();
  for (auto const &idx : indices)
    mean_r += r[idx];
  mean_r /= indices.size();

  real dr = 0.0;
  for (auto const &idx : indices)
    dr += norm_squared(r[idx] - mean_r);
  return sqrt(dr / indices.size());
}

void make_report::add_cur_scalars() const {
  auto &traj_div = rep->out_root.find<ioxx::sl4::div>(*traj_idx);
  auto &scalars = traj_div.find<ioxx::sl4::table>("scalars");
  auto &row = scalars->append_row();

  nitro::vector<int> all_indices(r.size());
  for (int idx = 0; idx < r.size(); ++idx)
    all_indices[idx] = idx;

  auto full_gyr = gyration(r, all_indices.get_view());
  auto full_Rg = gyration_radius(all_indices.get_view()),
       full_W = full_gyr.asphericity();
  auto full_rmsd = rmsd(all_indices.get_view());
  auto K = kinetic_energy(), E = *V + K;
  auto L = norm(r[r.size() - 1] - r[0]);

  row["TIME[tau]"] = quantity(*t).in("tau");
  row["EPOT[eps]"] = quantity(*V).in("eps");
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

void make_report::add_cur_model() const {
  auto xmd_model = *model;
  for (int idx = 0; idx < xmd_model.residues.size(); ++idx) {
    auto *ref_res = model->residues[idx].get();
    auto res_idx = res_map->at(ref_res);
    xmd_model.residues[idx]->pos = r[res_idx];
  }

  auto &pdb_model = rep->full_pdb.find_or_add_model(rep->model_serial);
  pdb_model = pdb_file(xmd_model).primary_model();
  pdb_model.model_serial = rep->model_serial;
  ++rep->model_serial;
}

void make_report::emit_pdb() const {
  auto pdb_path = with_ext(prefix, ".pdb");
  auto pdb_of = open_file(pdb_path);
  pdb_of << rep->full_pdb;
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
  auto &traj_div = rep->map_root.find<ioxx::sl4::div>(*traj_idx);
  auto &cur_div = traj_div.add<ioxx::sl4::div>();

  cur_div.add<ioxx::sl4::comment>("T = ", *t);

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

        auto dist = norm(r[cont.i1()] - r[cont.i2()]);
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

        auto dist = norm(r[cont.i1()] - r[cont.i2()]);
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

        auto dist = norm(r[bundle.i1()] - r[bundle.i2()]);
        pid_row["dist[A]"] = quantity(dist).in("A");
      }
    }

    num_comment = ioxx::sl4::comment("n = ", num_active_contacts);
  }

  cur_div.add<ioxx::sl4::comment>("angles");
  auto &angles_tab = cur_div.add<ioxx::sl4::table>();
  for (auto const &angle : model->angles) {
    auto &row = angles_tab->append_row();

    auto i1 = res_map->at(angle.res1), i2 = res_map->at(angle.res2),
         i3 = res_map->at(angle.res3);
    row["i1"] = i1;
    row["i2"] = i2;
    row["i3"] = i3;

    auto r1 = r[i1], r2 = r[i2], r3 = r[i3];
    row["theta"] = angle_value(r1, r2, r3);

    if (angle.theta.has_value())
      row["native"] = angle.theta.value();
  }

  cur_div.add<ioxx::sl4::comment>("dihedrals");
  auto &dih_tab = cur_div.add<ioxx::sl4::table>();
  for (auto const &dih : model->dihedrals) {
    auto &row = dih_tab->append_row();

    auto i1 = res_map->at(dih.res1), i2 = res_map->at(dih.res2),
         i3 = res_map->at(dih.res3), i4 = res_map->at(dih.res4);
    row["i1"] = i1;
    row["i2"] = i2;
    row["i3"] = i3;
    row["i4"] = i4;

    auto r1 = r[i1], r2 = r[i2], r3 = r[i3], r4 = r[i4];
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