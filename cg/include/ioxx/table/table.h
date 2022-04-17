#pragma once
#include <ioxx/convert.h>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace ioxx::table {
class table;
class columns;

class cell {
public:
  cell() = default;
  cell(cell const &other) = default;
  cell &operator=(cell const &other) = default;

  template <typename T> cell(T const &init) { *this = init; }

  template <typename T> cell &operator=(T const &init) {
    value_repr = convert<std::string>(init);
    return *this;
  }

  template <typename T> T as() const { return convert<T>(value_repr); }

public:
  std::string value_repr;
};

class row {
public:
  row(table const &super);
  row(table const &super, std::vector<std::string> const &values);

  row &operator=(std::vector<std::string> const &values);

  template <typename... Values>
  row(table const &super, Values &&...values) : super{super} {
    cells = {cell(values)...};
  }

  cell &operator[](std::string const &colname);
  cell const &operator[](std::string const &colname) const;

public:
  table const &super;
  std::vector<cell> cells;
};

class columns {
public:
  columns(table &super);
  columns(table &super, std::vector<std::string> const &names);

  columns &operator=(std::vector<std::string> const &names);

  void add(std::string const &name, std::string const &init = "");

  template <typename T> void add(std::string const &name, T const &init) {
    add(name, convert<std::string>(init));
  }

  int operator[](std::string const &name) const;

  int size() const;

public:
  table &super;
  std::vector<std::string> names;
  std::unordered_map<std::string, int> name_to_idx;
};

class table {
public:
  table();
  table(std::vector<std::string> const &colnames);
  table(std::initializer_list<std::string> const &colnames);

  row &add();
  row &add_(std::vector<std::string> const &values);

  template <typename... Values> row &add(Values &&...values) {
    auto values_list =
        std::vector<std::string>{convert<std::string>(values)...};
    return add_(values_list);
  }

  int num_cols() const;

public:
  std::vector<row> rows;
  columns cols;
};
} // namespace ioxx::table