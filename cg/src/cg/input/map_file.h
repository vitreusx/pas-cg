#pragma once
#include <filesystem>
#include <ioxx/csv.h>
#include <ioxx/ioxx.h>
#include <string>
#include <vector>

namespace cg {
class map_file {
public:
  struct contact {
    int i1, i2;
    double length;
    void connect(ioxx::row_proxy &proxy);
  };
  std::vector<contact> contacts;

  struct angle {
    int i1, i2, i3;
    double theta;
    void connect(ioxx::row_proxy &proxy);
  };
  std::vector<angle> angles;

  struct dihedral {
    int i1, i2, i3, i4;
    double phi;
    void connect(ioxx::row_proxy &proxy);
  };
  std::vector<dihedral> dihedrals;

  void shift(int shift_val);

  void link(ioxx::xyaml::proxy &proxy);
};
} // namespace cg