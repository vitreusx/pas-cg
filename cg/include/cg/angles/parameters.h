#pragma once
#include "aa_heur_pair.h"
#include <cg/files/files.h>
#include <cg/utils/quantity.h>
#include <optional>

namespace cg::angles {
struct heur_ang_parameters {
  bool enabled;
  struct pair_coeffs {
    char type1, type2;
    quantity poly[7];
  };
  std::unordered_map<aa_heur_pair, pair_coeffs> coeffs;
  void load(ioxx::xyaml::node const &node);
};

struct heur_dih_parameters {
  bool enabled;
  struct pair_coeffs {
    char type2, type3;
    quantity const_, sin, cos, sin2, cos2, sin_cos;
  };
  std::unordered_map<aa_heur_pair, pair_coeffs> coeffs;
  void load(ioxx::xyaml::node const &node);
};

struct nat_ang_parameters {
  bool enabled;
  quantity CBA;
  void load(ioxx::xyaml::node const &node);
};

struct nat_dih_parameters {
  bool enabled;
  std::string variant;

  struct complex_variant_parameters {
    quantity CDA, CDB;
    void load(ioxx::xyaml::node const &node);
  };
  complex_variant_parameters complex;

  struct simple_variant_parameters {
    quantity CDH;
    void load(ioxx::xyaml::node const &node);
  };
  simple_variant_parameters simple;

  void load(ioxx::xyaml::node const &node);
};

struct parameters {
  bool all_enabled, dih_enabled;
  heur_ang_parameters heur_ang;
  nat_ang_parameters nat_ang;
  heur_dih_parameters heur_dih;
  nat_dih_parameters nat_dih;
  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::angles