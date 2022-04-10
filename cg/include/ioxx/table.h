#pragma once
#include <ioxx/convert.h>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace ioxx::table {
class cell_proxy;
class row;
class row_proxy;
class table;

enum class proxy_mode { READING, WRITING };

class cell {
public:
  template <typename T> cell &operator=(T const &x) {
    value = convert<std::string>(x);
  }

  template <typename T> cell &operator<<(T const &x) { *this = x; }

  template <typename T> T as() const { return convert<T>(value); }

  template <typename T> cell &operator>>(T &x) const { x = this->as<T>(); }

public:
  std::string value;
};

class cell_proxy {
public:
  explicit cell_proxy(cell &cell_ref);
  explicit cell_proxy(cell const &cell_ref);

  template <typename T> cell_proxy &operator&(T &x) {
    switch (mode) {
    case proxy_mode::READING:
      cell_ref >> x;
      break;
    case proxy_mode::WRITING:
      cell_ref << x;
      break;
    }
    return *this;
  }

public:
  cell &cell_ref;
  proxy_mode mode;
};

class columns {
public:
  std::string const &operator[](size_t idx) const;
  size_t operator[](std::string const &col) const;

  std::string &operator[](size_t idx);
  size_t operator[](std::string const &col);

  struct iter {
    iter(int cur, int max);

    int cur, max;
    iter operator++(int);
    iter &operator++();
    int operator*() const;
    bool operator!=(iter const &other);
  };

  size_t count() const;
  iter begin() const;
  iter end() const;

public:
  std::vector<std::string> idx_to_name;
  std::unordered_map<std::string, int> name_to_idx;
};

class row {
public:
  row(columns &cols);

  cell &operator[](size_t idx);
  cell const &operator[](size_t idx) const;

  cell &operator[](std::string const &col);
  cell const &operator[](std::string const &col) const;

  template <typename T> row &operator=(T const &value);
  template <typename T> row &operator<<(T const &value);

  template <typename T> T as(T init = T()) const;
  template <typename T> row const &operator>>(T &value) const;

public:
  std::vector<cell> cells;
  columns &cols;
};

class row_proxy {
public:
  explicit row_proxy(row &row_ref);
  explicit row_proxy(row const &row_ref);

  cell_proxy operator[](size_t idx) const;
  cell_proxy operator[](std::string const &col) const;

public:
  row &row_ref;
  proxy_mode mode;
};

template <typename T> struct user_conv;

template <typename T> struct conv {
  void map(row_proxy &proxy, T &value) { user_conv<T>().map(proxy, value); }

  void map(row_proxy &proxy, T const &value) {
    map(proxy, const_cast<T &>(value));
  }
};

template <typename T> row &row::operator=(T const &value) {
  conv<T>().map(*this, value);
  return *this;
}

template <typename T> row &row::operator<<(T const &value) {
  conv<T>().map(*this, value);
  return *this;
}

template <typename T> T row::as(T init) const {
  conv<T>().map(*this, init);
  return std::move(init);
}

template <typename T> row const &row::operator>>(T &value) const {
  conv<T>().map(*this, value);
  return *this;
}

class table_csv_repr;

std::istream &operator>>(std::istream &is, table_csv_repr &repr);
std::ostream &operator<<(std::ostream &os, table_csv_repr const &repr);

class table_csv_repr {
public:
  friend std::istream &operator>>(std::istream &is, table_csv_repr &repr);
  friend std::ostream &operator<<(std::ostream &os, table_csv_repr const &repr);

public:
  table &table;
};

class table_sl4_repr;

std::istream &operator>>(std::istream &is, table_sl4_repr &repr);
std::ostream &operator<<(std::ostream &os, table_sl4_repr const &repr);

class table_sl4_repr {
public:
  friend std::istream &operator>>(std::istream &is, table_sl4_repr &repr);
  friend std::ostream &operator<<(std::ostream &os, table_sl4_repr const &repr);

public:
  table &table;
};

class table {
public:
  row &operator[](size_t idx);
  cell &operator[](std::pair<size_t, size_t> row_col_idx);
  cell &operator[](std::pair<size_t, std::string> row_name);
  row &append();

  row const &operator[](size_t idx) const;
  cell const &operator[](std::pair<size_t, size_t> row_col_idx) const;
  cell const &operator[](std::pair<size_t, std::string> row_name) const;

  table_csv_repr csv();
  table_sl4_repr sl4();

  template <typename T> struct row_iter {
    row_iter(std::vector<row> const &rows, int cur, int max)
        : rows{rows}, cur{cur}, max{max} {}

    std::vector<row> const &rows;
    int cur, max;

    row_iter operator++(int) {
      auto copy = *this;
      cur++;
      return copy;
    }

    row_iter &operator++() {
      ++cur;
      return *this;
    }

    T operator*() const { return rows[cur].as<T>(); }

    bool operator!=(row_iter const &other) { return cur != other.cur; }
  };

  template <typename T> struct rows_seq {
    rows_seq(std::vector<row> const &rows) : rows{rows} {}

    std::vector<row> const &rows;

    row_iter<T> begin() { return row_iter<T>(rows, 0, rows.size()); }

    row_iter<T> end() { return row_iter<T>(rows, rows.size(), rows.size()); }
  };

  template <typename T> rows_seq<T> view_as() const {
    return rows_seq<T>(rows);
  }

public:
  columns cols;
  std::vector<row> rows;
};

} // namespace ioxx::table