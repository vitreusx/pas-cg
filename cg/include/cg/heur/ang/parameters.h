#pragma once
#include <cg/heur/aa_heur_pair.h>
#include <cg/utils/quantity.h>
#include <ioxx/csv.h>
#include <ioxx/ioxx.h>
#include <unordered_map>

namespace cg::heur_ang {
struct parameters {
  bool enabled;

  struct pair_coeffs {
    char type1, type2;
    quantity poly[7];
    void connect(ioxx::row_proxy &proxy);
  };
  std::unordered_map<aa_heur_pair, pair_coeffs> coeffs;

  void load(ioxx::xyaml::node const &node);
};
} // namespace cg::heur_ang