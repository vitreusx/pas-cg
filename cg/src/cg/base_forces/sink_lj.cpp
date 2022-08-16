#include <cg/base_forces/per_pair.h>
#include <cg/base_forces/sink_lj.h>

namespace cg {
void sink_lj_specs::load(const ioxx::xyaml::node &node) {
  node["depth"] >> depth;
  node["r_low"] >> r_low;
  node["r_high"] >> r_high;
}

sink_lj_specs::operator sink_lj() const {
  return sink_lj(depth.value(), r_low.value(), r_high.value());
}
sink_lj_specs::sink_lj_specs(const lj_specs &specs)
    : depth(specs.depth), r_low(specs.r_min), r_high(specs.r_min) {}

void ss_sink_lj_specs::load(const ioxx::xyaml::node &node) {
  for (auto const &aa1 : amino_acid::all()) {
    for (auto const &aa2 : amino_acid::all()) {
      ss_specs[{aa1, aa2}] = {};
    }
  }

  if (auto def_n = node["default"]; def_n) {
    for (auto const &aa1 : amino_acid::all()) {
      for (auto const &aa2 : amino_acid::all()) {
        def_n["depth"] >> ss_specs[{aa1, aa2}].depth;
        def_n["r_low"] >> ss_specs[{aa1, aa2}].r_low;
        def_n["r_high"] >> ss_specs[{aa1, aa2}].r_high;
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

    if (auto r_low_n = per_pair_n["r_low"]; r_low_n) {
      auto r_low_csv = r_low_n.as<per_pair_csv<quantity>>();
      for (auto const &aa1 : amino_acid::all())
        for (auto const &aa2 : amino_acid::all())
          ss_specs[{aa1, aa2}].r_low = r_low_csv[{aa1, aa2}].assumed("A");
    }

    if (auto r_high_n = per_pair_n["r_high"]; r_high_n) {
      auto r_high_csv = r_high_n.as<per_pair_csv<quantity>>();
      for (auto const &aa1 : amino_acid::all())
        for (auto const &aa2 : amino_acid::all())
          ss_specs[{aa1, aa2}].r_high = r_high_csv[{aa1, aa2}].assumed("A");
    }
  }
}
} // namespace cg