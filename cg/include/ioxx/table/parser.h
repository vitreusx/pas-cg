#pragma once
#include "table.h"

namespace ioxx::table {
class parser {
public:
  virtual ~parser() {}
  virtual table read(std::istream &is) const = 0;
  virtual void write(std::ostream &os, table const &tab) const = 0;
};
} // namespace ioxx::table