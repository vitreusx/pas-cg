#pragma once
#include "parser.h"

namespace ioxx::table {
class sl4_parser : public parser {
public:
  virtual ~sl4_parser() = default;
  virtual table read(std::istream &is) const override;
  virtual void write(std::ostream &os, table const &tab) const override;

  void fit(table const &tab);
  void write(std::ostream &os, columns const &cols, bool comment = true) const;
  void write(std::ostream &os, row const &row, bool comment = false) const;

public:
  bool header = true;
  std::vector<dtype> dtypes;
  std::vector<int> col_width;
};
} // namespace ioxx::table