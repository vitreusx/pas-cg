#include <algorithm>
#include <cfloat>
#include <cg/files/table/sl4.h>
#include <cg/utils/text.h>
#include <cinttypes>

#define STR2(s) #s
#define STR(s) STR2(s)

namespace ioxx::table {
table sl4_parser::read(std::istream &) const {
  throw std::runtime_error("reading in SL4 format is not implemented");
}

static std::string sl4_format(std::string const &repr, dtype dt) {
  if (repr.empty())
    return "";

  switch (dt.t) {
  case dtype::tag::bool_: {
    auto value = convert<bool>(repr);
    return value ? "T" : "F";
  }
  case dtype::tag::uint_: {
    auto value = convert<uint64_t>(repr);
    return cg::format("%" PRIu64, value);
  }
  case dtype::tag::int_: {
    auto value = convert<int64_t>(repr);
    return cg::format("%" PRId64, value);
  }
  case dtype::tag::float_: {
    auto value = convert<double>(repr);
    return cg::format("%g", value);
  }
  case dtype::tag::string_:
  default: {
    return repr;
  }
  }
}

static void write_row(std::ostream &os, bool comment,
                      std::vector<std::string> const &fields) {
  if (comment)
    os << '#';
  else
    os << ' ';

  bool first = true;
  for (auto const &field : fields) {
    if (!first)
      os << "   ";
    os << field;
    first = false;
  }
}

static std::string uppercase(std::string s) {
  for (auto &c : s)
    c = std::toupper(c);
  return s;
}

void sl4_parser::write(std::ostream &os, const table &tab) const {
  auto copy_ = *this;
  copy_.fit(tab);

  copy_.write(os, tab.cols);
  for (auto const &row : tab.rows)
    copy_.write(os, row);
}

void sl4_parser::fit(const table &tab) {
  dtypes = tab.cols.dtypes;
  col_width = std::vector<int>(tab.num_cols(), 0);

  for (int col_idx = 0; col_idx < tab.num_cols(); ++col_idx) {
    auto dt = tab.cols.dtypes[col_idx];

    int width = 0;
    width = std::max(width, (int)tab.cols.names[col_idx].size());
    for (int row_idx = 0; row_idx < (int)tab.rows.size(); ++row_idx) {
      auto const &repr = tab.rows[row_idx].fields[col_idx];
      auto out_repr = sl4_format(repr, dt);
      width = std::max(width, (int)out_repr.size());
    }
    col_width[col_idx] = width;
  }
}

void sl4_parser::write(std::ostream &os, const columns &cols,
                       bool comment) const {
  std::vector<std::string> out_header(cols.size());

  for (int col_idx = 0; col_idx < (int)cols.size(); ++col_idx) {
    auto name = cols.names[col_idx];
    name = uppercase(name);
    auto padding = std::string(col_width[col_idx] - name.size(), ' ');
    name = padding + name;
    out_header[col_idx] = name;
  }

  write_row(os, comment, out_header);
}

void sl4_parser::write(std::ostream &os, row const &row, bool comment) const {
  std::vector<std::string> out_fields(col_width.size());

  for (int col_idx = 0; col_idx < (int)col_width.size(); ++col_idx) {
    auto dt = dtypes[col_idx];
    auto field = row.fields[col_idx];
    field = sl4_format(field, dt);
    auto padding = std::string(col_width[col_idx] - field.size(), ' ');
    field = padding + field;
    out_fields[col_idx] = field;
  }

  write_row(os, comment, out_fields);
}
} // namespace ioxx::table