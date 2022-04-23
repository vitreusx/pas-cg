#include <cg/files/table.h>

namespace ioxx::table {
dtype::dtype(tag t) : t{t} {}

cell_ref::cell_ref(table &super, int row_idx, int col_idx)
    : super{super}, row_idx{row_idx}, col_idx{col_idx} {}

std::string &cell_ref::field_ref() {
  return super.rows[row_idx].fields[col_idx];
}

std::string const &cell_ref::field_ref() const {
  return super.rows[row_idx].fields[col_idx];
}

dtype &cell_ref::dt() { return super.cols.dtypes[col_idx]; }

dtype const &cell_ref::dt() const { return super.cols.dtypes[col_idx]; }

cell_ref &cell_ref::operator=(cell_ref const &other) {
  field_ref() = super.rows[other.row_idx].fields[other.col_idx];
  dt() = other.dt();
  return *this;
}

cell_ref &cell_ref::operator=(cell_ref &&other) {
  field_ref() = other.field_ref();
  dt() = other.dt();
  return *this;
}

cell_const_ref::cell_const_ref(table const &super, int row_idx, int col_idx)
    : super{super}, row_idx{row_idx}, col_idx{col_idx} {}

cell_const_ref::cell_const_ref(const cell_ref &ref)
    : super{ref.super}, row_idx{ref.row_idx}, col_idx{ref.col_idx} {}

std::string const &cell_const_ref::field_ref() const {
  return super.rows[row_idx].fields[col_idx];
}

dtype const &cell_const_ref::dt() const { return super.cols.dtypes[col_idx]; }

row::row(table &super, int row_idx) : super{super}, row_idx{row_idx} {
  fields = std::vector<std::string>(super.num_cols());
}

row::row(table &super, int row_idx, const std::vector<std::string> &values)
    : row{super, row_idx} {
  *this = values;
}

row &row::operator=(const row &other) {
  row_idx = other.row_idx;
  fields = other.fields;
  return *this;
}

row &row::operator=(const std::vector<std::string> &values) {
  for (int col_idx = 0; col_idx < super.num_cols(); ++col_idx)
    (*this)[super.cols.names[col_idx]] = values[col_idx];
  return *this;
}

cell_ref row::operator[](const std::string &colname) {
  auto col_idx = super.cols[colname];
  return cell_ref(super, row_idx, col_idx);
}

cell_const_ref row::operator[](const std::string &colname) const {
  auto col_idx = static_cast<columns const &>(super.cols)[colname];
  return cell_const_ref(super, row_idx, col_idx);
}

columns::columns(table &super) : super{super} {}

columns::columns(table &super, const std::vector<std::string> &colnames)
    : super{super} {
  *this = colnames;
}

columns &columns::operator=(const columns &other) {
  names = other.names;
  dtypes = other.dtypes;
  name_to_idx = other.name_to_idx;
  return *this;
}

columns &columns::operator=(const std::vector<std::string> &colnames) {
  for (auto &row : super.rows)
    row = {};

  for (auto const &name : colnames)
    add(name);

  return *this;
}

int columns::operator[](const std::string &name) {
  if (name_to_idx.find(name) == name_to_idx.end())
    add(name);
  return name_to_idx.at(name);
}

int columns::operator[](const std::string &name) const {
  if (name_to_idx.find(name) == name_to_idx.end())
    throw std::runtime_error("column with name \"" + name + "\" not found");
  return name_to_idx.at(name);
}

int columns::size() const { return names.size(); }

void columns::add(const std::string &name) {
  auto col_idx = names.size();
  names.push_back(name);
  dtypes.emplace_back();
  name_to_idx[name] = col_idx;

  for (auto &row : super.rows)
    row.fields.emplace_back("");
}

columns::iter_t columns::begin() const { return names.begin(); }

columns::iter_t columns::end() const { return names.end(); }

table::table() : cols{*this} {}

table::table(const std::vector<std::string> &colnames)
    : cols{*this, colnames} {}

table::table(const std::initializer_list<std::string> &colnames)
    : cols{*this, colnames} {}

table &table::operator=(const table &other) {
  cols = other.cols;
  rows = {};
  for (auto const &row : other.rows)
    append_row() = row;
  return *this;
}

row &table::append_row() {
  int row_idx = rows.size();
  return rows.emplace_back(*this, row_idx);
}

int table::num_rows() const { return rows.size(); }

int table::num_cols() const { return cols.size(); }

} // namespace ioxx::table