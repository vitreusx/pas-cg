#pragma once
#include <ioxx/convert.h>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace ioxx {
class raw_csv_row;
template <typename Row = raw_csv_row> class csv;
class row_proxy;
class csv_header;

template <typename T> struct csv_connection {
  void operator()(row_proxy &proxy, T &value) const { value.connect(proxy); }
};

std::ostream &operator<<(std::ostream &os, csv_header const &header);

class csv_header {
public:
  csv_header() = default;
  explicit csv_header(row_proxy const &proxy);
  csv_header(std::initializer_list<std::string> col_names);
  explicit csv_header(std::vector<std::string> const& col_names);

  void insert(size_t idx, std::string const &name);

  size_t size() const;
  std::string const &operator[](size_t idx) const;
  size_t operator[](std::string const &s) const;

  friend std::ostream &operator<<(std::ostream &os, csv_header const &header);

  std::vector<std::string> const &col_names() const;

private:
  std::vector<std::string> idx_to_name;
  std::unordered_map<std::string, size_t> name_to_idx;
};

class raw_csv_cell {
public:
  explicit raw_csv_cell(std::string &ref);
  explicit raw_csv_cell(std::string const &ref);

  raw_csv_cell(raw_csv_cell const &other) = default;
  raw_csv_cell &operator=(raw_csv_cell const &other);

  template <typename T> raw_csv_cell &operator=(T const &value) {
    to_ref() = convert<std::string, T>(value);
    return *this;
  }

  template <typename T> raw_csv_cell& operator<<(T const& value) {
    return (*this = value);
  }

  raw_csv_cell(raw_csv_cell &&other) noexcept = default;
  raw_csv_cell &operator=(raw_csv_cell &&other) noexcept;

  template <typename T> T as() {
    return convert<T, std::string>(to_const_ref());
  }

private:
  std::string &to_ref();
  std::string const &to_const_ref() const;

  std::variant<std::string *, std::string const *> ptr;
};

class raw_csv_row {
public:
  raw_csv_row() = default;
  explicit raw_csv_row(csv_header const *header);
  explicit raw_csv_row(csv_header const *header, std::string const &line);

  raw_csv_cell operator[](size_t idx);
  raw_csv_cell operator[](size_t idx) const;
  raw_csv_cell operator[](std::string const &name);
  raw_csv_cell operator[](std::string const &name) const;

  csv_header const *header;
  std::vector<std::string> values;

  void connect(row_proxy &row);
};

enum class row_proxy_mode { LOAD, SAVE };

class cell_proxy : public raw_csv_cell {
public:
  explicit cell_proxy(raw_csv_cell const &base, row_proxy_mode mode);

  template <typename T> cell_proxy &operator<<(T const &value) {
    this->raw_csv_cell::operator=(value);
    return *this;
  }

  template <typename T> cell_proxy &operator>>(T &value) {
    value = this->as<T>();
    return *this;
  }

  template <typename T> cell_proxy &operator&(T &value) {
    switch (mode) {
    case row_proxy_mode::LOAD:
      *this >> value;
      break;
    case row_proxy_mode::SAVE:
      *this << value;
      break;
    }
    return *this;
  }

private:
  row_proxy_mode mode;
};

std::ostream &operator<<(std::ostream &os, row_proxy const &row);

class row_proxy : public raw_csv_row {
public:
  row_proxy(raw_csv_row const &data, row_proxy_mode mode);

  cell_proxy operator[](size_t idx);
  cell_proxy operator[](std::string const &name);

  friend std::ostream &operator<<(std::ostream &os, row_proxy const &row);

private:
  row_proxy_mode mode;

  friend class raw_csv_row;
  friend class csv_header;
};

template <typename Row> class csv {
public:
  csv() = default;

  explicit csv(std::string const &source, bool load_header = true) {
    std::stringstream ss{};
    ss << source;
    load(ss, load_header);
  }

  explicit csv(std::istream &file, bool load_header = true) {
    load(file, load_header);
  }

  template<typename U = Row, typename std::enable_if_t<std::is_same_v<U, raw_csv_row>, int> = 0>
  raw_csv_row& emplace_back() {
    return rows.emplace_back(&header.value());
  }

  std::string save(bool save_header = true) const {
    std::stringstream ss{};
    save(ss, save_header);
    return ss.str();
  }

  std::ostream &save(std::ostream &os, bool save_header = true) const {
    bool first = true;

    if (save_header && header.has_value()) {
      if (!first)
        os << '\n';
      os << header.value();
      first = false;
    }

    for (auto const &schema : rows) {
      row_proxy row(raw_csv_row(header_ptr()), row_proxy_mode::SAVE);
      csv_connection<Row>()(row, const_cast<Row &>(schema));

      if (!first)
        os << '\n';
      os << row;
      first = false;
    }

    return os;
  }

  std::optional<csv_header> header;
  std::vector<Row> rows;

  friend std::istream &operator>>(std::istream &is, csv<Row> const &csv_file) {
    csv_file.load(is, true);
    return is;
  }

  friend std::ostream &operator<<(std::ostream &os, csv<Row> const &csv_file) {
    csv_file.save(os, true);
    return os;
  }

private:
  csv_header const *header_ptr() const {
    if (header.has_value())
      return &header.value();
    else
      return nullptr;
  }

  void load(std::istream &is, bool load_header) {
    rows.clear();
    for (std::string line; std::getline(is, line);) {
      auto row_data = raw_csv_row(header_ptr(), line);
      auto proxy = row_proxy(row_data, row_proxy_mode::LOAD);

      if (load_header) {
        header = csv_header(proxy);
        load_header = false;
      } else {
        auto &row = rows.emplace_back();
        csv_connection<Row>()(proxy, row);
      }
    }
  }
};

} // namespace ioxx
