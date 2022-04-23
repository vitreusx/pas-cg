#include <cg/files/xyaml/file.h>
#include <cg/files/xyaml/table.h>

namespace ioxx::xyaml {
table_file::table_file(const std::filesystem::path &path,
                       const table::table &tab)
    : path{path}, tab{tab} {}

table::table *table_file::operator->() { return &tab; }

void table_file::load(const node &from) {
  file f;
  from >> f;

  path = f.rel_path;
  tab = table::csv_parser().load(f.fetch());
}

void table_file::save(node &to) const {
  file f;
  f.rel_path = path;
  f.source = table::csv_parser().dump(tab);
  to << f;
}

void user_repr<ioxx::table::table>::load(const node &from,
                                         ioxx::table::table &to) const {
  table_file f = from.as<table_file>();
  to = f.tab;
}

void user_repr<ioxx::table::table>::save(node &to,
                                         const ioxx::table::table &from) const {
  table_file f;
  f.tab = from;
  to << f;
}
} // namespace ioxx::xyaml