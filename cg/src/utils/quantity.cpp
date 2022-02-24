#include "utils/quantity.h"
#include <algorithm>
#include <cparse/builtin-features.inc>
#include <cparse/shunting-yard.h>
#include <stdexcept>
#include <unordered_map>
using namespace cg;
using namespace ioxx;

struct cparse_unit_init {
  cparse_unit_init() {
    using namespace cparse;
    cparse_startup();

    TokenMap &g = TokenMap::default_global();

    std::unordered_map<std::string, double> unit_map = {
        {"f77unit", f77unit}, {"A", angstrom},
        {"nm", nanometer},    {"m", meter},
        {"ns", nanosecond},   {"tau", tau},
        {"1/tau", 1 / tau},   {"micros", microsecond},
        {"ms", millisecond},  {"s", second},
        {"atom", atom},       {"mol", mol},
        {"eps", eps},         {"kcal", kcal},
        {"J", Joule},         {"kB", kB},
        {"K", Kelvin},        {"kg", kg},
        {"amu", amu},         {"f77mass", f77mass},
        {"e", echarge},       {"C", Coulomb},
        {"Amp", Ampere},      {"c", cspeed},
        {"H", Henry},         {"mu_0", mu_0},
        {"eps_0", eps_0},     {"rad", rad},
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

quantity::quantity() {
  this->mag = this->unit_value = 1.0;
  this->unit_str = "";
}

quantity::quantity(double value) {
  this->mag = value;
  this->unit_str = "";
  this->unit_value = 1.0;
}

quantity::quantity(double mag, std::string const &unit) {
  this->mag = mag;
  this->unit_str = unit;
  this->unit_value = parse(unit);
}

quantity::quantity(const std::string &repr) {
  auto space = std::find(repr.begin(), repr.end(), ' ');
  if (space == repr.end()) {
    mag = parse(repr);
    unit_str = "";
    unit_value = 1.0;
  } else {
    auto mag_str = std::string(repr.begin(), space);
    mag = parse(mag_str);
    unit_str = std::string(space + 1, repr.end());
    unit_value = parse(unit_str);
  }
}

quantity::operator double() const { return mag * unit_value; }

std::string quantity::repr() const {
  std::stringstream mag_ss{};
  mag_ss << std::scientific << mag;
  auto mag_repr = mag_ss.str();

  if (unit_str.empty())
    return mag_repr;
  else
    return mag_repr + " " + unit_str;
}

quantity &quantity::operator=(const quantity &other) {
  mag = other.mag;
  if (!other.unit_str.empty()) {
    unit_str = other.unit_str;
    unit_value = other.unit_value;
  }
  return *this;
}

quantity &quantity::operator=(quantity &&other) noexcept {
  *this = static_cast<quantity const &>(other);
  return *this;
}

quantity quantity::in(std::string const &new_unit) const {
  auto copy = *this;
  return copy.in_(new_unit);
}

quantity &quantity::in_(std::string const &new_unit) {
  double value = *this;
  unit_str = new_unit;
  if (!unit_str.empty())
    unit_value = parse(unit_str);
  else
    unit_value = 1.0;
  mag = value / unit_value;
  return *this;
}

std::istream &cg::operator>>(std::istream &is, quantity &value) {
  std::string repr;
  is >> repr;
  value = quantity(repr);
  return is;
}

std::ostream &cg::operator<<(std::ostream &os, quantity const &value) {
  os << value.repr();
  return os;
}
