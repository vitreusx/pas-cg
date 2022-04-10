#include <ioxx/table.h>

namespace ioxx::table {
cell_proxy::cell_proxy(cell &cell_ref)
    : cell_ref{cell_ref}, mode{proxy_mode::WRITING} {}

cell_proxy::cell_proxy(const cell &cell_ref)
    : cell_ref{const_cast<cell &>(cell_ref)}, mode{proxy_mode::READING} {}

std::string const &columns::operator[](size_t idx) const {
  if (idx >= idx_to_name.size())
    throw std::runtime_error("column index out of bounds");
  return idx_to_name.at(idx);
}

size_t columns::operator[](std::string const &col) const {
  if (name_to_idx.find(col) == name_to_idx.end())
    throw std::runtime_error("column not found");
  return name_to_idx.at(col);
}

std::string &columns::operator[](size_t idx) {
  if (idx >= count())
    idx_to_name.resize(idx + 1);
  return idx_to_name[idx];
}

size_t columns::operator[](std::string const &col) {
  if (name_to_idx.find(col) == name_to_idx.end()) {
    auto idx = count();
    name_to_idx[col] = idx;
    idx_to_name.push_back(col);
  }
  return name_to_idx.at(col);
}

columns::iter::iter(int cur, int max) : cur{cur}, max{max} {}

columns::iter &columns::iter::operator++() {
  ++cur;
  return *this;
}

columns::iter columns::iter::operator++(int) {
  auto copy = *this;
  cur++;
  return copy;
}

int columns::iter::operator*() const { return cur; }

bool columns::iter::operator!=(const iter &other) { return cur != other.cur; }

size_t columns::count() const { return idx_to_name.size(); }

columns::iter columns::begin() const { return columns::iter(0, count()); }

columns::iter columns::end() const { return columns::iter(count(), count()); }

row::row(columns &cols) : cols{cols} {}

cell &row::operator[](size_t idx) {
  if (idx >= cells.size())
    cells.resize(idx + 1);
  return cells[idx];
}

cell const &row::operator[](size_t idx) const {
  if (idx >= cells.size())
    throw std::runtime_error("column index out of bounds");
  return cells[idx];
}

cell &row::operator[](const std::string &col) {
  int idx = cols[col];
  return (*this)[idx];
}

cell const &row::operator[](const std::string &col) const {
  int idx = static_cast<columns const &>(cols)[col];
  return (*this)[idx];
}

row_proxy::row_proxy(row &row_ref)
    : row_ref{row_ref}, mode{proxy_mode::WRITING} {}

row_proxy::row_proxy(row const &row_ref)
    : row_ref{const_cast<row &>(row_ref)}, mode{proxy_mode::READING} {}

cell_proxy row_proxy::operator[](size_t idx) const {
  switch (mode) {
  case proxy_mode::READING:
    return cell_proxy(static_cast<row const &>(row_ref)[idx]);
  case proxy_mode::WRITING:
    return cell_proxy(row_ref[idx]);
  }
}

cell_proxy row_proxy::operator[](const std::string &col) const {
  switch (mode) {
  case proxy_mode::READING:
    return cell_proxy(static_cast<row const &>(row_ref)[col]);
  case proxy_mode::WRITING:
    return cell_proxy(row_ref[col]);
  }
}

row &table::operator[](size_t idx) {
  if (idx >= rows.size())
    throw std::runtime_error("row index out of bounds");
  return rows[idx];
}

cell &table::operator[](std::pair<size_t, size_t> row_col_idx) {
  auto const &[row_idx, col_idx] = row_col_idx;
  return (*this)[row_idx][col_idx];
}

cell &table::operator[](std::pair<size_t, std::string> row_name) {
  auto const &[row_idx, col_name] = row_name;
  return (*this)[row_idx][col_name];
}

row &table::append() { return rows.emplace_back(cols); }

row const &table::operator[](size_t idx) const {
  if (idx >= rows.size())
    throw std::runtime_error("row index out of bounds");
  return rows[idx];
}

cell const &table::operator[](std::pair<size_t, size_t> row_col_idx) const {
  auto const &[row_idx, col_idx] = row_col_idx;
  return (*this)[row_idx][col_idx];
}

cell const &table::operator[](std::pair<size_t, std::string> row_name) const {
  auto const &[row_idx, col_name] = row_name;
  return (*this)[row_idx][col_name];
}

} // namespace ioxx::table