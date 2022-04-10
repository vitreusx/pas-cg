#include <cg/base_forces/lj_variants.h>
#include <cg/utils/ioxx_interop.h>
#include <cg/utils/quantity.h>
namespace cg {

struct per_pair_csv {
  amino_acid type;
  std::unordered_map<amino_acid, quantity> values;

  void connect(ioxx::row_proxy &row) {
    type = amino_acid(row["type"].as<std::string>());
    for (auto const &col : row.header->col_names()) {
      if (col != "type")
        row[col] & values[amino_acid(col)];
    }
  }
};

void lj_variants::load(ioxx::xyaml::node const &node) {
  bb.r_min() = node["bb"]["r_min"].as<quantity>();
  bb.depth() = node["bb"]["depth"].as<quantity>();

  bs.r_min() = node["bs"]["r_min"].as<quantity>();
  bs.depth() = node["bs"]["depth"].as<quantity>();
  sb = bs;

  std::optional<quantity> def_depth, def_r_min, def_r_max;
  if (auto ss_def_values = node["ss"]["default"]; ss_def_values) {
    ss_def_values["depth"] >> def_depth;
    ss_def_values["r_min"] >> def_r_min;
    ss_def_values["r_max"] >> def_r_max;
  }

  std::optional<ioxx::csv<per_pair_csv>> per_pair_depth, per_pair_r_min,
      per_pair_r_max;
  if (auto ss_per_pair = node["ss"]["per pair"]; ss_per_pair) {
    ss_per_pair["depth"] >> per_pair_depth;
    ss_per_pair["r_min"] >> per_pair_r_min;
    ss_per_pair["r_max"] >> per_pair_r_max;
  }

  auto use_sinking = node["ss"]["use sinking variant"].as<bool>();
  if (use_sinking) {
    if (!def_r_min.has_value())
      def_r_min = bb.r_min();
  } else {
    def_r_max = def_r_min;
    if (per_pair_r_min.has_value())
      per_pair_r_max = per_pair_r_min;
    else if (per_pair_r_max.has_value())
      per_pair_r_min = per_pair_r_max;
  }

  for (auto const &aa1 : amino_acid::all()) {
    for (auto const &aa2 : amino_acid::all()) {
      ss[{aa1, aa2}] = sink_lj(bb);

      if (def_r_min.has_value())
        ss[{aa1, aa2}].r_min() = def_r_min.value();
      if (def_r_max.has_value())
        ss[{aa1, aa2}].r_max() = def_r_max.value();
      if (def_depth.has_value())
        ss[{aa1, aa2}].depth() = def_depth.value();
    }
  }

  if (per_pair_r_min.has_value()) {
    for (auto const &row : per_pair_r_min.value().rows) {
      auto aa1 = row.type;
      for (auto const &[aa2, r_min] : row.values) {
        ss[{aa1, aa2}].r_min() = r_min.assumed("A");
      }
    }
  }

  if (per_pair_r_max.has_value()) {
    for (auto const &row : per_pair_r_max.value().rows) {
      auto aa1 = row.type;
      for (auto const &[aa2, r_max] : row.values) {
        ss[{aa1, aa2}].r_max() = r_max.assumed("A");
      }
    }
  }

  if (per_pair_depth.has_value()) {
    for (auto const &row : per_pair_depth.value().rows) {
      auto aa1 = row.type;
      for (auto const &[aa2, depth] : row.values) {
        ss[{aa1, aa2}].depth() = depth.assumed("eps");
      }
    }
  }
}
} // namespace cg