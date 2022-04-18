#include <cg/utils/text.h>
#include <inttypes.h>
#include <ioxx/table/sl4.h>

namespace ioxx::table {
table sl4_parser::read(std::istream &is) const {
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
  case dtype::tag::string_: {
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
      os << ' ';
    os << field;
    first = false;
  }
}

void sl4_parser::write(std::ostream &os, const table &tab) const {
  std::vector<int> col_width(tab.num_cols(), 0);
  std::vector<std::string> out_header(tab.num_cols());
  std::vector<std::vector<std::string>> out_fields(
      tab.rows.size(), std::vector<std::string>(tab.num_cols()));

  for (int col_idx = 0; col_idx < tab.num_cols(); ++col_idx) {
    auto dt = tab.cols.dtypes[col_idx];
    out_header[col_idx] = tab.cols.names[col_idx];

    int width = 0;
    width = std::max(width, (int)tab.cols.names[col_idx].size());
    for (int row_idx = 0; row_idx < tab.rows.size(); ++row_idx) {
      auto const &repr = tab.rows[row_idx].fields[col_idx];
      auto out_repr = sl4_format(repr, dt);
      out_fields[row_idx][col_idx] = out_repr;
      width = std::max(width, (int)out_repr.size());
    }
    col_width[col_idx] = width;
  }

  for (int col_idx = 0; col_idx < tab.num_cols(); ++col_idx) {
    auto &out_header_col = out_header[col_idx];
    auto padding =
        std::string(col_width[col_idx] - (int)out_header_col.size(), ' ');
    out_header_col = padding + out_header_col;
  }

  for (int row_idx = 0; row_idx < tab.rows.size(); ++row_idx) {
    for (int col_idx = 0; col_idx < tab.num_cols(); ++col_idx) {
      auto &out_repr = out_fields[row_idx][col_idx];
      auto padding =
          std::string(col_width[col_idx] - (int)out_repr.size(), ' ');
      out_repr = padding + out_repr;
    }
  }

  write_row(os, true, out_header);
  for (auto const &row : out_fields) {
    os << '\n';
    write_row(os, false, row);
  }
}
} // namespace ioxx::table