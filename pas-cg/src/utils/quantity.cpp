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

    auto u = [](auto x) -> auto { return x.asDouble(); };

    /* Distance */
    auto f77unit = u(g["f77unit"] = 1.0);
    auto angstrom = u(g["angstrom"] = f77unit / 5.0);
    auto nanometer = u(g["nm"] = 10.0 * angstrom);
    auto meter = u(g["m"] = nanometer * 1.0e9);

    /* Time */
    auto nanosecond = u(g["ns"] = 1.0);
    auto tau = u(g["tau"] = 1.0 * nanosecond);
    auto microsecond = u(g["micros"] = nanosecond * 1.0e3);
    auto millisecond = u(g["ms"] = nanosecond * 1.0e6);
    auto second = u(g["s"] = nanosecond * 1.0e9);
    auto tau_inv = u(g["1/tau"] = 1.0 / tau);

    /* Quantity */
    auto atom = u(g["atom"] = 1.0);
    auto mol = u(g["mol"] = 6.02214076e23 / atom);

    /* Energy */
    auto eps = u(g["eps"] = 1.0); /* \approx 1.5kcal/mol */
    auto kcal = u(g["kcal"] = eps * mol / 1.57824959);
    auto Joule = u(g["J"] = kcal / 4184.0);

    /* Temperature */
    auto eps_kB = u(g["eps/kB"] = 1.0);
    auto kB = u(g["kB"] = eps / eps_kB);
    auto Kelvin = u(g["K"] = 1.380649e-23 * Joule / kB);

    /* Mass */
    auto kg = u(g["kg"] = Joule * second * second / (meter * meter));
    auto amu = u(g["amu"] = kg * 0.99999999965e-3 / mol);

    /**
     * In the Fortran version of the model, distance of \p f77unit, time of
     * \p tau, energy of \p eps and the average mass of an aminoacid were units;
     * these are however incongruent, it's a confirmed bug.
     */
    auto f77mass = u(g["f77mass"] = eps * tau * tau / (f77unit * f77unit));

    /* EM stuff */
    auto echarge = u(g["e"] = 1.0);
    auto Coulomb = u(g["C"] = echarge / 1.602176634e-19);
    auto Ampere = u(g["Amp"] = Coulomb / second);
    auto cspeed = u(g["c"] = 299792458.0 * meter / second);
    auto Henry =
        u(g["H"] = kg * meter * meter / (second * second * Ampere * Ampere));
    auto mu_0 = u(g["mu_0"] = 1.25663706212e-6 * Henry / meter);
    auto eps_0 = u(g["eps_0"] = 1.0 / (mu_0 * cspeed * cspeed));

    /* Degrees */
    auto rad = u(g["rad"] = 1.0);
    auto deg = u(g["deg"] = (2.0 * M_PI / 360.0) * rad);
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

quantity &quantity::in_(std::string const &new_unit) {
  double value = *this;
  unit_str = new_unit;
  unit_value = parse(unit_str);
  mag = value / unit_value;
  return *this;
}

std::string convert_impl<std::string, quantity>::impl(const quantity &q) {
  return q.repr();
}

quantity convert_impl<quantity, std::string>::impl(const std::string &s) {
  return quantity(s);
}