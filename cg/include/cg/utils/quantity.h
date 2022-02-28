#pragma once
#include "units.h"
#include <cmath>
#include <ioxx/convert.h>
#include <optional>
#include <string>

namespace cg {

class quantity {
public:
  quantity() = default;
  quantity(double internal_value);
  explicit quantity(double numerical_value, std::string const &unit);
  explicit quantity(std::string const &repr);

  quantity(quantity const &other) = default;
  quantity(quantity &&other) = default;

  quantity &operator=(quantity const &other);
  quantity &operator=(quantity &&other) noexcept;

  friend std::istream &operator>>(std::istream &is, quantity &value);
  friend std::ostream &operator<<(std::ostream &os, quantity const &value);

  quantity &assumed_(std::string const &def_unit);
  quantity assumed(std::string const &def_unit) const;
  double in(std::string const &unit) const;

  operator double() const;
  double num_value() const;
  std::string repr() const;

private:
  double numerical_value = 1.0;
  std::optional<std::string> unit_str = std::nullopt;
  std::optional<std::string> _def_unit = std::nullopt;
};
} // namespace cg
