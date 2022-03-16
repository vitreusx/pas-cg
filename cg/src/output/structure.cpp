#include "output/structure.h"
#include "utils/quantity.h"
#include <ioxx/ioxx.h>
namespace cg::out {

static auto tether_value(vec3r r1, vec3r r2) { return norm(r2 - r1); }

struct tether_row {
  int i1, i2;
  double dist;
  std::optional<double> nat_dist;

  void connect(ioxx::row_proxy &row) {
    row["i1"] << i1;
    row["i2"] << i2;
    row["dist[A]"] << quantity(dist).in("A");
    if (nat_dist.has_value())
      row["nat_dist[A]"] << quantity(nat_dist.value()).in("A");
  }
};

static auto angle_value(cg::vec3r r1, cg::vec3r r2, cg::vec3r r3) {
  auto r12_u = unit(r2 - r1), r23_u = unit(r3 - r2);
  return acos(dot(r12_u, r23_u));
}

struct angle_row {
  int i1, i2, i3;
  double theta;
  std::optional<double> nat_theta;

  void connect(ioxx::row_proxy &row) {
    row["i1"] << i1;
    row["i2"] << i2;
    row["i3"] << i3;
    row["theta[rad]"] << quantity(theta).in("rad");
    if (nat_theta.has_value())
      row["nat_theta[rad]"] << quantity(nat_theta.value()).in("rad");
  }
};

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

struct dihedral_row {
  int i1, i2, i3, i4;
  double phi;
  std::optional<double> nat_phi;

  void connect(ioxx::row_proxy &row) {
    row["i1"] << i1;
    row["i2"] << i2;
    row["i3"] << i3;
    row["i4"] << i4;
    row["phi[rad]"] << quantity(phi).in("rad");
    if (nat_phi.has_value())
      row["nat_phi[rad]"] << quantity(nat_phi.value()).in("rad");
  }
};

void add_structure::report_to(report_data &report) const {
  if (report.report_files) {
    ioxx::xyaml::csv<tether_row> tether_file;
    tether_file.path = "tethers.csv";
    tether_file.data.header = {"i1", "i2", "dist[A]", "nat_dist[A]"};
    for (auto const &tether : model->tethers) {
      tether_row row;
      row.i1 = res_map->at(tether.res1);
      row.i2 = res_map->at(tether.res2);
      row.dist = tether_value(r[row.i1], r[row.i2]);
      row.nat_dist = tether.length;
      tether_file.data.rows.push_back(row);
    }
    report.for_snap["tethers"] = tether_file;

    ioxx::xyaml::csv<angle_row> angle_file;
    angle_file.path = "angles.csv";
    angle_file.data.header = {"i1", "i2", "i3", "theta[rad]", "nat_theta[rad]"};
    for (auto const &angle : model->angles) {
      angle_row row;
      row.i1 = res_map->at(angle.res1);
      row.i2 = res_map->at(angle.res2);
      row.i3 = res_map->at(angle.res3);
      row.theta = angle_value(r[row.i1], r[row.i2], r[row.i3]);
      row.nat_theta = angle.theta;
      angle_file.data.rows.push_back(row);
    }
    report.for_snap["angles"] = angle_file;

    ioxx::xyaml::csv<dihedral_row> dih_file;
    dih_file.path = "dihedrals.csv";
    dih_file.data.header = {"i1", "i2", "i3", "i4", "phi[rad]", "nat_phi[rad]"};
    for (auto const &dih : model->dihedrals) {
      dihedral_row row;
      row.i1 = res_map->at(dih.res1);
      row.i2 = res_map->at(dih.res2);
      row.i3 = res_map->at(dih.res3);
      row.i4 = res_map->at(dih.res4);
      row.phi = dihedral_value(r[row.i1], r[row.i2], r[row.i3], r[row.i4]);
      row.nat_phi = dih.phi;
      dih_file.data.rows.push_back(row);
    }
    report.for_snap["dihedrals"] = dih_file;
  }
}
} // namespace cg::out