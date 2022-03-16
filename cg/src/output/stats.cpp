#include "output/stats.h"
#include "utils/ioxx_interop.h"
#include "utils/math.h"
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <ioxx/ioxx.h>

namespace cg::out {

quantity gyration_stuff::radius_of_gyration() const { return norm(lambda); }

quantity gyration_stuff::asphericity() const {
  return (real)1.5 * pow(lambda.x(), 2.0) - (real)0.5 * norm_squared(lambda);
}

gyration_stuff add_stats::gyration(const nitro::vector<int> &indices) const {
  Eigen::MatrixX3d positions(indices.size(), 3);
  for (int idx = 0; idx < indices.size(); ++idx)
    positions.row(idx) = convert<double>(r[indices[idx]]).transpose();

  Eigen::RowVector3d center = positions.colwise().mean();
  positions.rowwise() -= center;

  Eigen::Vector3d lambda =
      Eigen::JacobiSVD<decltype(positions)>(positions).singularValues();
  lambda = lambda.cwiseSqrt();

  gyration_stuff data;
  data.lambda = vec3r(lambda);
  return data;
}

real add_stats::rmsd(const nitro::vector<int> &indices) const {
  auto n = indices.size();
  Eigen::MatrixX3d orig_R(n, 3), cur_R(n, 3);
  for (int idx = 0; idx < n; ++idx) {
    orig_R.row(idx) = convert<double>(orig_r[indices[idx]]).transpose();
    cur_R.row(idx) = convert<double>(r[indices[idx]]).transpose();
  }

  orig_R -= orig_R.colwise().mean();
  cur_R -= cur_R.colwise().mean();

  Eigen::Matrix3d H = orig_R.transpose() * cur_R;
  auto svd = Eigen::JacobiSVD<decltype(H)>(H);
  auto svdU = svd.matrixU(), svdV = svd.matrixV();
  auto d = sign((svdV.transpose() * svdU).determinant());
  auto R = svdV * Eigen::Vector3d(1, 1, d).asDiagonal() * svdU.transpose();

  auto rmsd = sqrt((cur_R * R - orig_R).rowwise().squaredNorm().sum());
  return (real)rmsd;
}

void add_stats::report_to(report_data &report) const {
  auto stats_node = report.for_snap["stats"];
  stats_node["t"] = *t;
  report.add_step_scalar("t", *t);
  stats_node["V"] = *V;
  report.add_step_scalar("V", *V);

  real K = 0;
  for (int idx = 0; idx < v.size(); ++idx)
    K += (real)0.5 * mass[(uint8_t)atype[idx]] * ipow<2>(norm(v[idx]));

  report.add_step_scalar("K", stats_node["K"], K);
  report.add_step_scalar("E", stats_node["E"], K + *V);

  auto all_indices = nitro::vector<int>(r.size());
  for (int idx = 0; idx < r.size(); ++idx)
    all_indices[idx] = idx;

  auto gyr_total = gyration(all_indices);
  report.add_step_scalar("Rg", stats_node["Rg"],
                         gyr_total.radius_of_gyration().in("A"));
  report.add_step_scalar("W", stats_node["W"],
                         gyr_total.asphericity().in("A**2"));

  if (report.report_files) {
    ioxx::xyaml::csv<> chains_csv;
    chains_csv.path = "chains.csv";
    chains_csv.data.header = {"idx", "first", "last", "L", "Rg", "W"};

    int num_chains = chain_first.size();
    for (int idx = 0; idx < num_chains; ++idx) {
      auto &row = chains_csv.data.emplace_back();
      row["idx"] = idx;
      auto first = chain_first[idx], last = chain_last[idx];
      row["first"] = first;
      row["last"] = last;
      row["L"] = quantity(norm(r[last] - r[first])).in("A");

      nitro::vector<int> chain_indices;
      for (int res_idx = 0; res_idx < r.size(); ++res_idx)
        if (chain_idx[res_idx] == idx)
          chain_indices.push_back(res_idx);

      auto gyr_chain = gyration(chain_indices);
      row["Rg"] = gyr_chain.radius_of_gyration().in("A");
      row["W"] = gyr_chain.asphericity().in("A**2");
    }
    report.for_snap["chains"] = chains_csv;
  }
}
} // namespace cg::out