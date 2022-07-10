#include "per_pair.h"
#include <cg/base_forces/lj.h>
#include <cg/utils/hash.h>
#include <cg/utils/ioxx_interop.h>

namespace cg {
void lj_specs::load(const ioxx::xyaml::node &node) {
  node["depth"] >> depth;
  node["r_min"] >> r_min;
}

lj_specs::operator lj() const {
  return lj(depth.value(), r_min.value());
}

void ss_lj_specs::load(const ioxx::xyaml::node &node) {
  for (auto const &aa1 : amino_acid::all()) {
    for (auto const &aa2 : amino_acid::all()) {
      ss_specs[{aa1, aa2}] = {};
    }
  }

  if (auto def_n = node["default"]; def_n) {
    for (auto const &aa1 : amino_acid::all()) {
      for (auto const &aa2 : amino_acid::all()) {
        def_n["depth"] >> ss_specs[{aa1, aa2}].depth;
        def_n["r_min"] >> ss_specs[{aa1, aa2}].r_min;
      }
    }
  }

  if (auto per_pair_n = node["per pair"]; per_pair_n) {
    if (auto depth_n = per_pair_n["depth"]; depth_n) {
      auto depth_csv = depth_n.as<per_pair_csv<quantity>>();
      for (auto const &aa1 : amino_acid::all())
        for (auto const &aa2 : amino_acid::all())
          ss_specs[{aa1, aa2}].depth = depth_csv[{aa1, aa2}].assumed("eps");
    }

    if (auto r_min_n = per_pair_n["r_min"]; r_min_n) {
      auto r_min_csv = r_min_n.as<per_pair_csv<quantity>>();
      for (auto const &aa1 : amino_acid::all())
        for (auto const &aa2 : amino_acid::all())
          ss_specs[{aa1, aa2}].r_min = r_min_csv[{aa1, aa2}].assumed("A");
    }
  }
}
} // namespace cg