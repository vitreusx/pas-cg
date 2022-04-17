#pragma once
#include "parser.h"

namespace ioxx::table {
class csv_parser : public parser {
public:
  virtual ~csv_parser() = default;
  virtual table read(std::istream &is) const override;
  virtual void write(std::ostream &os, table const &tab) const override;

public:
  bool header = true;
};
} // namespace ioxx::table