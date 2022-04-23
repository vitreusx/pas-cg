#pragma once
#include "table.h"

namespace ioxx::table {
class parser {
public:
  virtual ~parser() {}
  virtual table read(std::istream &is) const = 0;
  virtual void write(std::ostream &os, table const &tab) const = 0;

  table load(std::string const &source) const;
  std::string dump(table const &tab) const;
};
} // namespace ioxx::table