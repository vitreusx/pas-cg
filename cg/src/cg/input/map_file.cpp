#include <cg/input/map_file.h>
#include <fstream>
#include <sstream>

namespace cg {
void map_file::shift(int shift_val) {
  for (auto &cont : contacts) {
    cont.i1 += shift_val;
    cont.i2 += shift_val;
  }

  for (auto &ang : angles) {
    ang.i1 += shift_val;
    ang.i2 += shift_val;
    ang.i3 += shift_val;
  }

  for (auto &dih : dihedrals) {
    dih.i1 += shift_val;
    dih.i2 += shift_val;
    dih.i3 += shift_val;
    dih.i4 += shift_val;
  }
}

void map_file::load(const ioxx::xyaml::node &n) {
  ioxx::table::table contacts_csv, angles_csv, dihedrals_csv;

  n["contacts"] >> contacts_csv;
  for (auto const &row : contacts_csv.rows) {
    auto &cont = contacts.emplace_back();
    row["i1"] >> cont.i1;
    row["i2"] >> cont.i2;
    row["length"] >> cont.length;
  }

  n["angles"] >> angles_csv;
  for (auto const &row : angles_csv.rows) {
    auto &ang = angles.emplace_back();
    row["i1"] >> ang.i1;
    row["i2"] >> ang.i2;
    row["i3"] >> ang.i3;
    row["theta"] >> ang.theta;
  }

  n["dihedrals"] >> dihedrals_csv;
  for (auto const &row : dihedrals_csv.rows) {
    auto &dih = dihedrals.emplace_back();
    row["i1"] >> dih.i1;
    row["i2"] >> dih.i2;
    row["i3"] >> dih.i3;
    row["i4"] >> dih.i4;
    row["phi"] >> dih.phi;
  }
}
} // namespace cg
