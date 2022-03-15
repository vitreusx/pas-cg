#include "utils/quantity.h"
#include <algorithm>
#include <cparse/builtin-features.inc>
#include <cparse/shunting-yard.h>
#include <regex>
#include <stdexcept>
#include <unordered_map>

namespace cg {

struct cparse_unit_init {
  cparse_unit_init() {
    using namespace cparse;
    cparse_startup();

    TokenMap &g = TokenMap::default_global();

    std::unordered_map<std::string, double> unit_map = {{"f77unit", f77unit},
                                                        {"A", angstrom},
                                                        {"nm", nanometer},
                                                        {"m", meter},
                                                        {"ns", nanosecond},
                                                        {"tau", tau},
                                                        {"micros", microsecond},
                                                        {"ms", millisecond},
                                                        {"s", second},
                                                        {"atom", atom},
                                                        {"mol", mol},
                                                        {"eps", eps},
                                                        {"kcal", kcal},
                                                        {"J", Joule},
                                                        {"kB", kB},
                                                        {"K", Kelvin},
                                                        {"kg", kg},
                                                        {"amu", amu},
                                                        {"f77mass", f77mass},
                                                        {"e", echarge},
                                                        {"C", Coulomb},
                                                        {"Amp", Ampere},
                                                        {"c", cspeed},
                                                        {"H", Henry},
                                                        {"mu_0", mu_0},
                                                        {"eps_0", eps_0},
                                                        {"rad", rad},
                                                        {"deg", deg}};

    for (auto const &[name, val] : unit_map) {
      g[name] = val;
    }
  }
};

static double parse(std::string const &s) {
  static cparse_unit_init init;
  return calculator::calculate(s.c_str()).asDouble();
}

quantity::quantity(double value) { numerical_value = value; }

quantity::quantity(double numerical_value, std::string const &unit) {
  this->numerical_value = numerical_value;
  this->unit_str = unit;
}

static bool is_number(std::string const &str) {
  static std::regex number_re("-?[0-9]+([\\.][0-9]+)?");
  return std::regex_match(str, number_re);
}

quantity::quantity(const std::string &repr) {
  auto space = std::find(repr.begin(), repr.end(), ' ');
  if (space == repr.end()) {
    auto value = parse(repr);
    if (is_number(repr)) {
      numerical_value = value;
    } else {
      numerical_value = 1.0;
      unit_str = repr;
    }
  } else {
    auto num_val_str = std::string(repr.begin(), space);
    numerical_value = parse(num_val_str);
    unit_str = std::string(space + 1, repr.end());
  }
}

quantity::operator double() const {
  auto unit_value = 1.0;
  if (unit_str.has_value())
    unit_value = parse(unit_str.value());
  else if (_def_unit.has_value())
    unit_value = parse(_def_unit.value());
  return numerical_value * unit_value;
}

std::string quantity::repr() const {
  std::stringstream num_val_ss{};
  num_val_ss << std::scientific << numerical_value;
  auto num_val_repr = num_val_ss.str();

  if (!unit_str.has_value())
    return num_val_repr;
  else
    return num_val_repr + " " + unit_str.value();
}

quantity &quantity::operator=(const quantity &other) {
  numerical_value = other.numerical_value;
  unit_str = other.unit_str;
  return *this;
}

quantity &quantity::operator=(quantity &&other) noexcept {
  *this = static_cast<quantity const &>(other);
  return *this;
}

std::istream &operator>>(std::istream &is, quantity &value) {
  std::string repr;
  is >> repr;
  value = quantity(repr);
  return is;
}

std::ostream &operator<<(std::ostream &os, quantity const &value) {
  os << value.repr();
  return os;
}

quantity &quantity::assumed_(const std::string &def_unit) {
  _def_unit = def_unit;
  return *this;
}

quantity quantity::assumed(const std::string &def_unit) const {
  auto copy = *this;
  copy.assumed_(def_unit);
  return copy;
}

double quantity::in(const std::string &unit) const {
  return (double)*this / parse(unit);
}

double quantity::num_value() const { return numerical_value; }
} // namespace cg
