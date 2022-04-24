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

  auto &for_traj = traj[traj_idx];
  for_traj.traj_idx = traj_idx;
}

void make_report::operator()() const {
  if (*t >= rep->stats_last + stats_every) {
    add_cur_scalars();
    emit_out();
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

void make_report::add_cur_scalars() const {
  auto &scalars = rep->traj[*traj_idx].scalars;
  auto &row = scalars.append_row();

  nitro::vector<int> all_indices(r.size());
  for (int idx = 0; idx < r.size(); ++idx)
    all_indices[idx] = idx;

  auto full_gyr = gyration(r, all_indices.get_view());
  auto full_Rg = full_gyr.radius(), full_W = full_gyr.asphericity();
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
    int num_active_contacts = 0;
    for (auto const &cont : *nc->contacts) {
      if (nc->is_active(cont))
        ++num_active_contacts;
    }

    row["ICN"] = num_active_contacts;
  }
}

void make_report::emit_out() const {
  ioxx::sl4::div out_root;
  for (auto const &[traj_idx_, for_traj] : rep->traj) {
    auto &traj_div = out_root.add<ioxx::sl4::div>();
    traj_div.add<ioxx::sl4::comment>("trajectory ", traj_idx_ + 1);
    traj_div.add<ioxx::sl4::table>(for_traj.scalars);
  }

  auto out_path = with_ext(prefix, ".out");
  auto out_of = open_file(out_path);
  out_of << out_root;
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

void make_report::emit_map() const {
  auto map_path = with_ext(prefix, ".map");
  auto map_of = open_file(map_path);
}

} // namespace cg::out