#include <algorithm>
#include <ioxx/table/csv.h>
#include <regex>

namespace ioxx::table {
enum class csv_parser_state { UnquotedField, QuotedField, QuotedQuote };

// Based on: https://stackoverflow.com/a/30338543
static std::vector<std::string> read_csv_row(std::string const &row) {
  auto state = csv_parser_state::UnquotedField;
  std::vector<std::string> fields;

  for (auto c : row) {
    switch (state) {
    case csv_parser_state::UnquotedField:
      switch (c) {
      case ',': // end of field
        fields.emplace_back("");
        break;
      case '"':
        state = csv_parser_state::QuotedField;
        break;
      default:
        fields.back().push_back(c);
        break;
      }
      break;
    case csv_parser_state::QuotedField:
      switch (c) {
      case '"':
        state = csv_parser_state::QuotedQuote;
        break;
      default:
        fields.back().push_back(c);
        break;
      }
      break;
    case csv_parser_state::QuotedQuote:
      switch (c) {
      case ',': // , after closing quote
        fields.emplace_back("");
        state = csv_parser_state::UnquotedField;
        break;
      case '"': // "" -> "
        fields.back().push_back('"');
        state = csv_parser_state::QuotedField;
        break;
      default: // end of quote
        state = csv_parser_state::UnquotedField;
        break;
      }
      break;
    }
  }

  return fields;
}

table csv_parser::read(std::istream &is) const {
  table res;

  bool read_header = false;
  for (std::string row; std::getline(is, row);) {
    auto fields = read_csv_row(row);
    if (header && !read_header) {
      res.cols = fields;
      read_header = true;
    } else {
      res.add_(fields);
    }
  }

  return res;
}

std::string escape(std::string const &value) {
  auto has_comma = (std::find(value.begin(), value.end(), ',') != value.end());
  auto has_quote = (std::find(value.begin(), value.end(), '"') != value.end());

  if (has_comma || has_quote) {
    static const std::regex quote_re("\"");
    std::string quoted_value;
    std::regex_replace(std::back_inserter(quoted_value), value.begin(),
                       value.end(), quote_re, "\"\"\"");
    quoted_value = "\"" + quoted_value + "\"";
    return quoted_value;
  } else {
    return value;
  }
}

template <typename T>
static void write_row(std::ostream &os, std::vector<T> const &values) {
  for (int cell_idx = 0; cell_idx < values.size(); ++cell_idx) {
    if (cell_idx > 0)
      os << ',';
    os << values[cell_idx];
  }
}

void csv_parser::write(std::ostream &os, const table &tab) const {
  if (header) {
    write_row(os, tab.cols.names);
    os << '\n';
  }

  for (int row_idx = 0; row_idx < tab.rows.size(); ++row_idx) {
    if (row_idx > 0)
      os << '\n';
    write_row(os, tab.rows[row_idx].fields);
  }
}
} // namespace ioxx::table