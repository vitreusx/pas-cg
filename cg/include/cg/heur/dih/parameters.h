#pragma once
#include <cg/heur/aa_heur_pair.h>
#include <cg/utils/quantity.h>
#include <ioxx/csv.h>
#include <ioxx/xyaml.h>
#include <unordered_map>

namespace cg::heur_dih {
struct parameters {
  bool enabled;

  struct pair_coeffs {
    char type2, type3;
    quantity const_, sin, cos, sin2, cos2, sin_cos;
    void connect(ioxx::row_proxy &proxy);
  };
  std::unordered_map<aa_heur_pair, pair_coeffs> coeffs;

  void connect(ioxx::xyaml_proxy &proxy);
};
} // namespace cg::heur_dih