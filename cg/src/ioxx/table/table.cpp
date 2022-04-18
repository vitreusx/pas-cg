#include <ioxx/table.h>

namespace ioxx::table {
dtype::dtype(tag t) : t{t} {}

cell_ref::cell_ref(std::string &field_ref, dtype &dt)
    : field_ref{field_ref}, dt{dt} {}

cell_ref &cell_ref::operator=(cell_ref const &other) {
  field_ref = other.field_ref;
  dt = other.dt;
  return *this;
}

cell_ref &cell_ref::operator=(cell_ref &&other) {
  field_ref = other.field_ref;
  dt = other.dt;
  return *this;
}

cell_const_ref::cell_const_ref(const std::string &field_ref)
    : field_ref{field_ref} {}

cell_const_ref::cell_const_ref(const cell_ref &ref)
    : field_ref{ref.field_ref} {}

row::row(table &super) : super{super} {
  fields = std::vector<std::string>(super.num_cols());
}

row::row(table &super, const std::vector<std::string> &values) : super{super} {
  *this = values;
}

row &row::operator=(const std::vector<std::string> &values) {
  fields = std::vector<std::string>(super.num_cols());
  for (int col_idx = 0; col_idx < super.num_cols(); ++col_idx)
    fields[col_idx] = values[col_idx];
  return *this;
}

cell_ref row::operator[](const std::string &colname) {
  auto col_idx = super.cols[colname];
  return cell_ref(fields[col_idx], super.cols.dtypes[col_idx]);
}

cell_const_ref row::operator[](const std::string &colname) const {
  return cell_const_ref(fields[super.cols[colname]]);
}

columns::columns(table &super) : super{super} {}

columns::columns(table &super, const std::vector<std::string> &names)
    : super{super} {
  *this = names;
}

columns &columns::operator=(const std::vector<std::string> &names) {
  this->names = names;
  dtypes = std::vector<dtype>(names.size());

  for (int col_idx = 0; col_idx < names.size(); ++col_idx)
    name_to_idx[names[col_idx]] = col_idx;

  std::vector<std::string> init(names.size());
  for (auto &row : super.rows)
    row = init;

  return *this;
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