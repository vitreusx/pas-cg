#include "output/rmsd.h"
#include "utils/ioxx_interop.h"
#include <Eigen/Dense>
#include <Eigen/SVD>

using namespace cg::out;

template <typename T> int sign(T const &value) {
  return (T(0) < value) - (value < T(0));
}

void compute_rmsd::report_to(report_state &report) const {
  int n = orig_r.size();
  Eigen::MatrixX3d orig_R(n, 3), cur_R(n, 3);
  for (int idx = 0; idx < n; ++idx) {
    orig_R.row(idx) = convert<double>(orig_r[idx]).transpose();
    cur_R.row(idx) = convert<double>(cur_r[idx]).transpose();
  }

  orig_R -= orig_R.colwise().mean();
  cur_R -= cur_R.colwise().mean();

  Eigen::Matrix3d H = orig_R.transpose() * cur_R;
  auto svd = Eigen::JacobiSVD<decltype(H)>(H);
  auto U = svd.matrixU(), V = svd.matrixV();
  auto d = sign((V.transpose() * U).determinant());
  auto R = V * Eigen::Vector3d(1, 1, d).asDiagonal() * U.transpose();

  auto rmsd = sqrt((cur_R * R - orig_R).rowwise().squaredNorm().sum());
  report.current["stats"]["root mean square deviation"] = rmsd;
}