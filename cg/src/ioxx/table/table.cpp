#include <ioxx/table.h>

namespace ioxx::table {
row::row(const table &super) : super{super} {
  cells = std::vector<cell>(super.num_cols());
}

row::row(const table &super, const std::vector<std::string> &values)
    : super{super} {
  *this = values;
}

row &row::operator=(const std::vector<std::string> &values) {
  cells = std::vector<cell>(super.num_cols());
  for (int col_idx = 0; col_idx < super.num_cols(); ++col_idx)
    cells[col_idx] = cell(values[col_idx]);
  return *this;
}

cell &row::operator[](const std::string &colname) {
  return cells[super.cols[colname]];
}

cell const &row::operator[](const std::string &colname) const {
  return cells[super.cols[colname]];
}

columns::columns(table &super) : super{super} {}

columns::columns(table &super, const std::vector<std::string> &names)
    : super{super} {
  *this = names;
}

columns &columns::operator=(const std::vector<std::string> &names) {
  this->names = names;
  for (int col_idx = 0; col_idx < names.size(); ++col_idx)
    name_to_idx[names[col_idx]] = col_idx;

  std::vector<std::string> init(names.size());
  for (auto &row : super.rows)
    row = init;

  return *this;
}

void columns::add(const std::string &name, const std::string &init) {
  auto col_idx = names.size();
  names.push_back(name);
  name_to_idx[name] = col_idx;

  for (auto &row : super.rows)
    row.cells.emplace_back(init);
}

int columns::operator[](const std::string &name) const {
  if (name_to_idx.find(name) == name_to_idx.end())
    throw std::runtime_error("column with name \"" + name + "\" not found");
  return name_to_idx.at(name);
}

int columns::size() const { return names.size(); }

table::table() : cols{*this} {}

table::table(const std::vector<std::string> &colnames)
    : cols{*this, colnames} {}

table::table(std::initializer_list<std::string> const &colnames)
    : cols{*this, colnames} {}

row &table::add() { return rows.emplace_back(*this); }

row &table::add_(const std::vector<std::string> &values) {
  return rows.emplace_back(*this, values);
}

int table::num_cols() const { return cols.size(); }

} // namespace ioxx::table