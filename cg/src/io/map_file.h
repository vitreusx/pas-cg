#pragma once
#include <filesystem>
#include <ioxx/csv.h>
#include <ioxx/xyaml.h>
#include <string>
#include <vector>

namespace cg {
class map_file {
public:
  struct contact {
    int i1, i2;
    double length;
    bool connect(ioxx::row_proxy &proxy);
  };
  std::vector<contact> contacts;

  struct angle {
    int i1, i2, i3;
    double theta;
    bool connect(ioxx::row_proxy &proxy);
  };
  std::vector<angle> angles;

  struct dihedral {
    int i1, i2, i3, i4;
    double phi;
    bool connect(ioxx::row_proxy &proxy);
  };
  std::vector<dihedral> dihedrals;

  void shift(int shift_val);

  bool connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg