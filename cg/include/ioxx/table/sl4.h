#pragma once
#include "parser.h"

namespace ioxx::table {
class sl4_parser : public parser {
public:
  virtual ~sl4_parser() = default;
  virtual table read(std::istream &is) const override;
  virtual void write(std::ostream &os, table const &tab) const override;

public:
  bool header = true;
};
} // namespace ioxx::table