#pragma once
#include <cmath>
#include <ioxx/convert.h>
#include <string>

namespace cg {

class quantity {
public:
  quantity();
  quantity(double value);
  explicit quantity(double mag, std::string const &unit);
  explicit quantity(std::string const &repr);

  quantity(quantity const &other) = default;
  quantity(quantity &&other) = default;

  quantity &operator=(quantity const &other);
  quantity &operator=(quantity &&other) noexcept;

  friend std::istream &operator>>(std::istream &is, quantity &value);
  friend std::ostream &operator<<(std::ostream &os, quantity const &value);

  operator double() const;
  std::string repr() const;

  quantity in(std::string const &new_unit) const;
  quantity &in_(std::string const &new_unit);

private:
  double mag, unit_value;
  std::string unit_str;
};
} // namespace cg
