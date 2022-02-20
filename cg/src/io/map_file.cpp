#include "map_file.h"
#include <fstream>
#include <ioxx/csv.h>
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

void map_file::contact::connect(ioxx::row_proxy &proxy) {
  proxy["i1"] & i1;
  proxy["i2"] & i2;
  proxy["length"] & length;
}

void map_file::angle::connect(ioxx::row_proxy &proxy) {
  proxy["i1"] & i1;
  proxy["i2"] & i2;
  proxy["i3"] & i3;
  proxy["theta"] & theta;
}

void map_file::dihedral::connect(ioxx::row_proxy &proxy) {
  proxy["i1"] & i1;
  proxy["i2"] & i2;
  proxy["i3"] & i3;
  proxy["i4"] & i4;
  proxy["phi"] & phi;
}

void map_file::connect(ioxx::xyaml_node_proxy &proxy) {
  ioxx::xyaml_csv_node<contact> contacts_csv;
  ioxx::xyaml_csv_node<angle> angles_csv;
  ioxx::xyaml_csv_node<dihedral> dihedrals_csv;

  if (!proxy.loading()) {
    contacts_csv.csv_file.rows = contacts;
    angles_csv.csv_file.rows = angles;
    dihedrals_csv.csv_file.rows = dihedrals;
  }

  proxy["contacts"] & contacts_csv;
  proxy["angles"] & angles_csv;
  proxy["dihedrals"] & dihedrals_csv;

  if (proxy.loading()) {
    contacts = contacts_csv.csv_file.rows;
    angles = angles_csv.csv_file.rows;
    dihedrals = dihedrals_csv.csv_file.rows;
  }
}
} // namespace cg
