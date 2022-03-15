#include "output/gyration.h"
#include "utils/ioxx_interop.h"
#include <Eigen/SVD>
namespace cg::out {

void report_gyration_stuff::report_to(report_data &report) const {
  Eigen::MatrixX3d positions(r.size(), 3);
  for (int idx = 0; idx < r.size(); ++idx)
    positions.row(idx) = convert<double>(r[idx]).transpose();

  Eigen::RowVector3d center = positions.colwise().mean();
  positions.rowwise() -= center;

  Eigen::Vector3d lambda =
      Eigen::JacobiSVD<decltype(positions)>(positions).singularValues();
  lambda = lambda.cwiseSqrt();

  auto node = report.for_step["kinematics"];
  node["lambda"] = vec3r(lambda);

  auto radius_of_gyration = quantity(lambda.norm()).in("A");
  report.add_step_scalar("Rg [A]", node["radius of gyration [A]"],
                         radius_of_gyration);
  auto asph =
      (real)1.5 * pow(lambda.x(), 2.0) - (real)0.5 * lambda.squaredNorm();
  auto asph_in_A = quantity(asph).in("A**2");
  report.add_step_scalar("W [A**2]", node["asphericity [A**2]"], asph_in_A);
}
} // namespace cg::out