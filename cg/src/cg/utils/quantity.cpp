#include <algorithm>
#include <cg/utils/quantity.h>
#include <cg/utils/units.h>
#include <cparse/builtin-features.inc>
#include <cparse/shunting-yard.h>
#include <mutex>
#include <regex>
#include <stdexcept>
#include <unordered_map>

namespace cg {

static packToken float_(TokenMap scope) {
  auto x = std::stof(scope["x"].asString());
  return (double)x;
}

class quantity_parser {
public:
  quantity_parser() {
    using namespace cparse;
    cparse_startup();

    unit_map = {{"f77unit", f77unit},
                {"A", angstrom},
                {"nm", nanometer},
                {"m", meter},
                {"ns", nanosecond},
                {"tau", tau},
                {"micros", microsecond},
                {"ms", millisecond},
                {"s", second},
                {"atom", atom},
                {"residue", atom},
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

    auto &token_map = TokenMap::default_global();
    token_map["float"] = CppFunction(&float_, {"x"}, "float");
  }

  double compute_p(std::string const &x) const {
    auto &token_map = TokenMap::default_global();
    for (auto const &[name, val] : unit_map)
      token_map[name] = val.p;
    return calculator::calculate(x.c_str()).asDouble();
  }

  double compute_q(std::string const &x) const {
    auto &token_map = TokenMap::default_global();
    for (auto const &[name, val] : unit_map)
      token_map[name] = val.q;
    return calculator::calculate(x.c_str()).asDouble();
  }

  ratio parse_unit(std::string const &x) const {
    auto p = compute_p(x), q = compute_q(x);
    return ratio(p, q);
  }

  double parse_num(std::string const &x) const {
    return calculator::calculate(x.c_str()).asDouble();
  }

private:
  std::unordered_map<std::string, ratio> unit_map;
};

static ratio parse_unit(std::string const &s) {
  static std::mutex mut;
  std::lock_guard<std::mutex> guard(mut);
  static quantity_parser parser;
  return parser.parse_unit(s);
}

static ratio parse_num(std::string const &s) {
  static std::mutex mut;
  std::lock_guard<std::mutex> guard(mut);
  static quantity_parser parser;
  return parser.parse_num(s);
}

quantity::quantity(ratio value) {
  value_ = value;
  unit_str = "";
}

quantity::quantity(ratio value, std::string const &unit) {
  value_ = value;
  unit_str = unit;
}

quantity::quantity(char const *repr) : quantity(std::string(repr)) {}

static bool is_number(std::string const &str) {
  static std::regex number_re("-?[0-9]+([\\.][0-9]+)?");
  return std::regex_match(str, number_re);
}

quantity::quantity(const std::string &repr) {
  auto space = std::find(repr.begin(), repr.end(), ' ');
  if (space == repr.end()) {
    if (is_number(repr)) {
      value_ = parse_num(repr);
      unit_str = "";
    } else {
      value_ = 1.0;
      unit_str = repr;
    }
  } else {
    auto num_val_str = std::string(repr.begin(), space);
    value_ = parse_num(num_val_str);
    unit_str = std::string(space + 1, repr.end());
  }
}

quantity::operator ratio() const {
  ratio unit_value = 1.0;
  if (!unit_str.empty())
    unit_value = parse_unit(unit_str);
  return value_ * unit_value;
}

std::string quantity::repr() const {
  std::stringstream num_val_ss{};
  num_val_ss << std::scientific << (double)value_;
  auto num_val_repr = num_val_ss.str();

  if (unit_str.empty())
    return num_val_repr;
  else
    return num_val_repr + " " + unit_str;
}

quantity::quantity(const quantity &other) {
  value_ = other.value_;
  if (def_unit_.has_value()) {
    unit_str = other.unit_str.empty() ? def_unit_.value() : other.unit_str;
  } else {
    unit_str = other.unit_str;
  }
}

quantity::quantity(quantity &&other)
    : quantity(static_cast<quantity const &>(other)) {}

quantity &quantity::operator=(const quantity &other) {
  value_ = other.value_;
  if (def_unit_.has_value()) {
    unit_str = other.unit_str.empty() ? def_unit_.value() : other.unit_str;
  } else {
    unit_str = other.unit_str;
  }
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

quantity &quantity::assumed(const std::string &def_unit) {
  def_unit_ = def_unit;
  return *this;
}

quantity quantity::assumed(const std::string &def_unit) const {
  return quantity(value_, unit_str.empty() ? def_unit : unit_str);
}

double quantity::value_in(const std::string &unit) const {
  return this->operator ratio() / parse_unit(unit);
}

ratio quantity::value() const {
  return value_;
}

quantity::operator double() const {
  return this->operator ratio();
}
} // namespace cg
