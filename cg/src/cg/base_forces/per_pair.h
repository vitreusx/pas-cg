#pragma once
#include <cg/amino/amino_acid.h>
#include <cg/files/files.h>
#include <cg/utils/hash.h>

namespace cg {
template <typename T> struct per_pair_csv {
  std::unordered_map<std::pair<amino_acid, amino_acid>, T> values;

  T const &operator[](std::pair<amino_acid, amino_acid> const &p) {
    return values.at(p);
  }

  void load(ioxx::xyaml::node const &n) {
    auto tab = n.as<ioxx::table::table>();
    for (auto const &row : tab.rows) {
      auto aa1 = amino_acid(row["type"].as<std::string>());
      for (auto const &col : tab.cols) {
        if (col != "type") {
          auto aa2 = amino_acid(col);
          auto val = row[col].as<T>();
          values[{aa1, aa2}] = val;
        }
      }
    }
  }
};
} // namespace cg