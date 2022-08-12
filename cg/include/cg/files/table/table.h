#pragma once
#include <cg/files/convert.h>
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
  enum class tag {
    uint_,
    int_,
    float_,
    bool_,
    string_
  };

  template <typename T> static dtype from_type() {
    if constexpr (std::is_same_v<T, bool>) {
      return dtype(tag::bool_);
    } else if constexpr (std::is_same_v<T, char> ||
                         std::is_same_v<T, unsigned char>) {
      return dtype(tag::string_);
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
  explicit cell_ref(table &super, int row_idx, int col_idx);

  cell_ref(cell_ref const &other) = delete;
  cell_ref(cell_ref &&other) = default;

  cell_ref &operator=(cell_ref const &other);
  cell_ref &operator=(cell_ref &&other);

  std::string &field_ref();
  std::string const &field_ref() const;

  dtype &dt();
  dtype const &dt() const;

  template <typename T> cell_ref &operator=(T const &value) {
    field_ref() = convert<std::string>(value);
    dt() = dtype::from_type<T>();
    return *this;
  }

  template <typename T> cell_ref &operator<<(T const &value) {
    return (*this = value);
  }

  template <typename T> T as() const {
    return convert<T>(field_ref());
  }

  template <typename T> cell_ref const &operator>>(T &value) const {
    value = as<T>();
    return *this;
  }

public:
  table &super;
  int row_idx, col_idx;
};

class cell_const_ref {
public:
  explicit cell_const_ref(table const &super, int row_idx, int col_idx);
  cell_const_ref(cell_ref const &ref);

  std::string const &field_ref() const;
  dtype const &dt() const;

  cell_const_ref(cell_const_ref const &other) = delete;
  cell_const_ref(cell_const_ref &&other) = default;

  cell_const_ref &operator=(cell_const_ref const &other) = delete;
  cell_const_ref &operator=(cell_const_ref &&other) = delete;

  template <typename T> T as() const {
    return convert<T>(field_ref());
  }

  template <typename T> cell_const_ref const &operator>>(T &value) const {
    value = as<T>();
    return *this;
  }

public:
  table const &super;
  int row_idx, col_idx;
};

class row {
public:
  row(table &super, int row_idx);
  row(table &super, int row_idx, std::vector<std::string> const &values);

  row(row const &other) = default;
  row &operator=(row const &other);

  row &operator=(std::vector<std::string> const &values);

  cell_ref operator[](std::string const &colname);
  cell_const_ref operator[](std::string const &colname) const;

public:
  table &super;
  int row_idx;
  std::vector<std::string> fields;
};

class columns {
public:
  explicit columns(table &super);
  explicit columns(table &super, std::vector<std::string> const &colnames);

  columns(columns const &other) = default;
  columns &operator=(columns const &other);

  columns &operator=(std::vector<std::string> const &colnames);

  void add(std::string const &name);

  int operator[](std::string const &name);
  int operator[](std::string const &name) const;

  int size() const;

public:
  table &super;
  std::vector<std::string> names;
  std::vector<dtype> dtypes;
  std::unordered_map<std::string, int> name_to_idx;

public:
  using iter_t = decltype(std::declval<decltype(names) const>().begin());
  iter_t begin() const;
  iter_t end() const;
};

class table {
public:
  table();
  table(std::vector<std::string> const &colnames);
  table(std::initializer_list<std::string> const &colnames);

  table(table const &other) = default;
  table &operator=(table const &other);

  row &append_row();
  int num_rows() const;
  int num_cols() const;

public:
  std::vector<row> rows;
  columns cols;
};
} // namespace ioxx::table