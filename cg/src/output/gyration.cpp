#include "output/gyration.h"
#include "utils/ioxx_interop.h"
#include <Eigen/SVD>
using namespace cg::out;

void report_gyration_stuff::report_to(report_state &report) const {
  Eigen::MatrixX3d positions(r.size(), 3);
  for (int idx = 0; idx < r.size(); ++idx)
    positions.row(idx) = convert<double>(r[idx]).transpose();

  Eigen::RowVector3d center_of_mass = positions.colwise().sum() / r.size();
  positions.rowwise() -= center_of_mass;

  Eigen::Vector3d lambda =
      Eigen::JacobiSVD<decltype(positions)>(positions).singularValues();

  auto node = report.current["kinematics"];
  node["lambda"] = vec3r(lambda);
  node["radius of gyration"] = lambda.norm();
  node["asphericity"] =
      (real)1.5 * pow(lambda.z(), 2.0) - (real)0.5 * lambda.squaredNorm();
}