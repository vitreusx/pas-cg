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

class dtype {
public:
  dtype() = default;
  enum class tag { uint_, int_, float_, bool_, string_ };

  template <typename T> static dtype from_type() {
    if constexpr (std::is_same_v<T, bool>) {
      return dtype(tag::bool_);
    } else if constexpr (std::is_integral_v<T>) {
      if constexpr (std::is_unsigned_v<T>)
        return dtype(tag::uint_);
      else
        return dtype(tag::int_);
    } else if constexpr (std::is_floating_point_v<T>) {
      return dtype(tag::float_);
    } else {
      return dtype(tag::string_);
    }
  }

public:
  explicit dtype(tag t);
  tag t = tag::string_;
};

class cell_ref {
public:
  explicit cell_ref(std::string &field_ref, dtype &dt);

  cell_ref(cell_ref const &other) = delete;
  cell_ref(cell_ref &&other) = default;

  cell_ref &operator=(cell_ref const &other);
  cell_ref &operator=(cell_ref &&other);

  template <typename T> cell_ref &operator=(T const &value) {
    field_ref = convert<std::string>(value);
    dt = dtype::from_type<T>();
    return *this;
  }

  template <typename T> T as() const { return convert<T>(field_ref); }

public:
  std::string &field_ref;
  dtype &dt;
};

class cell_const_ref {
public:
  explicit cell_const_ref(std::string const &field_ref);
  cell_const_ref(cell_ref const &ref);

  cell_const_ref(cell_const_ref const &other) = delete;
  cell_const_ref(cell_const_ref &&other) = default;

  cell_const_ref &operator=(cell_const_ref const &other) = delete;
  cell_const_ref &operator=(cell_const_ref &&other) = delete;

  template <typename T> T as() const { return convert<T>(field_ref); }

public:
  std::string const &field_ref;
};

class row {
public:
  row(table &super);
  row(table &super, std::vector<std::string> const &values);

  row &operator=(std::vector<std::string> const &values);

  template <typename... Values>
  row(table &super, Values &&...values) : super{super} {
    fields = {cell(values)...};
  }

  cell_ref operator[](std::string const &colname);
  cell_const_ref operator[](std::string const &colname) const;

public:
  table &super;
  std::vector<std::string> fields;
};

class columns {
public:
  columns(table &super);
  columns(table &super, std::vector<std::string> const &names);

  columns &operator=(std::vector<std::string> const &names);

  template <typename T> void add(std::string const &name, T const &init);

  int operator[](std::string const &name) const;

  int size() const;

public:
  table &super;
  std::vector<std::string> names;
  std::vector<dtype> dtypes;
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

template <typename T>
void columns::add(const std::string &name, const T &init) {
  auto col_idx = names.size();
  names.push_back(name);
  dtypes.push_back(dtype::from_type<T>());
  name_to_idx[name] = col_idx;

  for (auto &row : super.rows)
    row.fields.emplace_back(init);
}
} // namespace ioxx::table