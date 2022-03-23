#include "csv.h"
#include "read_csv_row.h"
using namespace ioxx;

csv_header::csv_header(std::initializer_list<std::string> col_names):
  csv_header(std::vector<std::string>(col_names)) {}

csv_header::csv_header(const std::vector<std::string> &col_names) {
  size_t col_idx = 0;
  for (auto const &col_name : col_names)
    insert(col_idx++, col_name);
}

csv_header::csv_header(const row_proxy &proxy) {
  size_t col_idx = 0;
  for (auto const &col_name : proxy.values)
    insert(col_idx++, col_name);
}

void csv_header::insert(size_t idx, std::string const &name) {
  if (idx + 1 > idx_to_name.size())
    idx_to_name.resize(idx + 1);

  idx_to_name[idx] = name;
  name_to_idx[name] = idx;
}

size_t csv_header::size() const { return idx_to_name.size(); }

std::string const &csv_header::operator[](size_t idx) const {
  return idx_to_name[idx];
}

size_t csv_header::operator[](std::string const &s) const {
  return name_to_idx.at(s);
}

std::ostream &ioxx::operator<<(std::ostream &os, csv_header const &header) {
  for (size_t idx = 0; idx < header.idx_to_name.size(); ++idx) {
    if (idx > 0)
      os << ',';
    os << header.idx_to_name[idx];
  }
  return os;
}
std::vector<std::string> const &csv_header::col_names() const {
  return idx_to_name;
}

cell_proxy::cell_proxy(raw_csv_cell const &base, row_proxy_mode mode)
    : raw_csv_cell{base}, mode{mode} {}

std::ostream &ioxx::operator<<(std::ostream &os, row_proxy const &row) {
  for (size_t idx = 0; idx < row.values.size(); ++idx) {
    if (idx > 0)
      os << ',';
    os << row.values[idx];
  }
  return os;
}

void raw_csv_row::connect(row_proxy &proxy) {
  switch (proxy.mode) {
  case row_proxy_mode::LOAD:
    *this = static_cast<raw_csv_row &>(proxy);
    break;
  case row_proxy_mode::SAVE:
    static_cast<raw_csv_row &>(proxy) = *this;
    break;
  }
}

raw_csv_cell raw_csv_row::operator[](size_t idx) {
  return raw_csv_cell(values[idx]);
}

raw_csv_cell raw_csv_row::operator[](size_t idx) const {
  return raw_csv_cell(values[idx]);
}

raw_csv_cell raw_csv_row::operator[](const std::string &name) {
  return (*this)[(*header)[name]];
}

raw_csv_cell raw_csv_row::operator[](const std::string &name) const {
  return (*this)[(*header)[name]];
}

raw_csv_row::raw_csv_row(const csv_header *header) : header{header} {
  if (header != nullptr)
    values = std::vector<std::string>(header->size());
}

raw_csv_row::raw_csv_row(const csv_header *header, const std::string &line)
    : header{header} {
  values = readCSVRow(line);
}

raw_csv_cell::raw_csv_cell(std::string &ref) : ptr{&ref} {}

raw_csv_cell::raw_csv_cell(const std::string &ref) : ptr{&ref} {}

raw_csv_cell &raw_csv_cell::operator=(const raw_csv_cell &other) {
  to_ref() = other.to_const_ref();
  return *this;
}

raw_csv_cell &raw_csv_cell::operator=(raw_csv_cell &&other) noexcept {
  to_ref() = std::move(other.to_ref());
  return *this;
}

std::string &raw_csv_cell::to_ref() { return *std::get<std::string *>(ptr); }

std::string const &raw_csv_cell::to_const_ref() const {
  if (std::holds_alternative<std::string *>(ptr))
    return *std::get<std::string *>(ptr);
  else
    return *std::get<std::string const *>(ptr);
}

row_proxy::row_proxy(const raw_csv_row &data, row_proxy_mode mode)
    : raw_csv_row(data), mode{mode} {}

cell_proxy row_proxy::operator[](size_t idx) {
  return cell_proxy(this->raw_csv_row::operator[](idx), mode);
}

cell_proxy row_proxy::operator[](std::string const &name) {
  return cell_proxy(this->raw_csv_row::operator[](name), mode);
}