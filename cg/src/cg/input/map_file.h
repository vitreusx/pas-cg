#pragma once
#include <cg/files/files.h>
#include <filesystem>
#include <string>
#include <vector>

namespace cg {
class map_file {
public:
  struct contact {
    int i1, i2;
    double length;
  };
  std::vector<contact> contacts;

  struct angle {
    int i1, i2, i3;
    double theta;
  };
  std::vector<angle> angles;

  struct dihedral {
    int i1, i2, i3, i4;
    double phi;
  };
  std::vector<dihedral> dihedrals;

  void shift(int shift_val);

  void load(ioxx::xyaml::node const &n);
};
} // namespace cg