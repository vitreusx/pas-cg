#pragma once
#include "../csv.h"
#include "file.h"

namespace ioxx::xyaml {
template <typename Row = raw_csv_row> class csv {
public:
  std::optional<std::filesystem::path> path;
  ioxx::csv<Row> data;
};

template <typename Row> struct user_repr<csv<Row>> {
  void load(node const &from, csv<Row> &to) const {
    auto csv_file = from.as<file>();
    to.path = csv_file.rel_path;
    to.data = ioxx::csv<Row>(csv_file.fetch());
  }

  void save(node &to, csv<Row> const &from) const {
    file csv_file;
    csv_file.rel_path = from.path;
    csv_file.source = from.data.save();
    to << csv_file;
  }
};

template <typename Row> struct user_repr<ioxx::csv<Row>> {
  void load(node const &from, ioxx::csv<Row> &to) const {
    auto csv_file = from.as<file>();
    to = ioxx::csv<Row>(csv_file.fetch());
  }

  void save(node &to, ioxx::csv<Row> const &from) const { to << from.save(); }
};

} // namespace ioxx::xyaml