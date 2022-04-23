#pragma once
#include "../table.h"
#include "node.h"

namespace ioxx::xyaml {
template <> struct user_repr<ioxx::table::table> {
  void load(node const &from, ioxx::table::table &to) const;
  void save(node &to, ioxx::table::table const &from) const;
};

class table_file {
public:
  table_file() = default;
  explicit table_file(std::filesystem::path const &path,
                      table::table const &tab);

  table::table *operator->();

  void load(node const &from);
  void save(node &to) const;

public:
  std::optional<std::filesystem::path> path;
  table::table tab;
};

} // namespace ioxx::xyaml