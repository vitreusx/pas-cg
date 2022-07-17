#pragma once

namespace cg {
class ratio {
public:
  constexpr ratio(double x) : ratio(x, 1.0) {}

  explicit constexpr ratio(double p, double q) : p{p}, q{q} {}

  friend constexpr ratio operator*(ratio const &x, ratio const &y) {
    return ratio(x.p * y.p, x.q * y.q);
  }

  friend constexpr ratio operator*(double x, ratio const &y) {
    return ratio(x * y.p, y.q);
  }

  friend constexpr ratio operator*(ratio const &x, double y) {
    return ratio(x.p * y, x.q);
  }

  constexpr ratio &operator*=(ratio const &y) {
    *this = *this * y;
    return *this;
  }

  friend constexpr ratio operator/(ratio const &x, ratio const &y) {
    return ratio(x.p * y.q, x.q * y.p);
  }

  friend constexpr ratio operator/(double x, ratio const &y) {
    return ratio(x * y.q, y.p);
  }

  friend constexpr ratio operator/(ratio const &x, double y) {
    return ratio(x.p, x.q * y);
  }

  constexpr ratio &operator/=(ratio const &y) {
    *this = *this / y;
    return *this;
  }

  constexpr operator double() const {
    return p / q;
  }

  double p, q;
};
} // namespace cg