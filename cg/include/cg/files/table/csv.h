#pragma once
#include "parser.h"

namespace ioxx::table {
class csv_parser : public parser {
public:
  ~csv_parser() override = default;
  table read(std::istream &is) const override;
  void write(std::ostream &os, table const &tab) const override;

public:
  bool header = true;
};
} // namespace ioxx::table