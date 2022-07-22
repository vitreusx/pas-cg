#pragma once
#include <cg/files/convert.h>
#include <cg/utils/ratio.h>
#include <cmath>
#include <optional>
#include <string>

namespace cg {
class quantity;

std::istream &operator>>(std::istream &is, quantity &value);
std::ostream &operator<<(std::ostream &os, quantity const &value);

class quantity {
public:
  quantity() = default;
  quantity(ratio value);
  explicit quantity(ratio value, std::string const &unit);
  explicit quantity(std::string const &repr);
  explicit quantity(char const *repr);

  template <typename T> quantity(T const &value) : quantity((ratio)value) {}

  template <typename T>
  quantity(T const &value, std::string const &unit)
      : quantity((ratio)value, unit) {}

  quantity(quantity const &other);
  quantity(quantity &&other);
  quantity &operator=(quantity const &other);
  quantity &operator=(quantity &&other) noexcept;

  friend std::istream &operator>>(std::istream &is, quantity &value);
  friend std::ostream &operator<<(std::ostream &os, quantity const &value);

  quantity &assumed(std::string const &def_unit);
  quantity assumed(std::string const &def_unit) const;

  double value_in(std::string const &unit) const;

  operator ratio() const;
  ratio value() const;
  std::string repr() const;

  inline operator double() const {
    return redux;
  }

private:
  ratio value_ = 1.0;
  double redux;
  std::string unit_str = "";
  std::optional<std::string> def_unit_;
};
} // namespace cg
