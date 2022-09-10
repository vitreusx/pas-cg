#pragma once
#include <cg/types/amp.h>
#include <cg/vect/vect.h>

namespace cg::pid {
enum lambda_version {
  COSINE,
  ALGEBRAIC
};

template <typename E>
struct lambda_expr {
  EXPR(psi_0, alpha, version)

  template <typename T>
  auto supp(T const &psi) const {
    return abs(alpha() * (psi - psi_0())) < (real)M_PI;
  }

  template <typename T>
  auto operator()(T const &psi) const {
    if constexpr (nitro::def::is_lane_like_v<T>) {
      return lane_call(psi);
    } else {
      return reg_call(psi);
    }
  }

private:
  template <typename T>
  auto lane_call(T const &psi) const {
    lambda_version ver;
    if constexpr (nitro::def::is_vcl_lane_v<std::decay_t<decltype(version())>>)
      ver = static_cast<lambda_version>(version()[0]);
    else
      ver = version();

    if (ver == COSINE) {
      auto s = alpha() * (psi - psi_0());
      auto mask = abs(s) < (real)M_PI;
      auto val = select(mask, (real)0.5f * cos(s) + (real)0.5f, (T)0);
      auto deriv = select(mask, -(real)0.5f * alpha() * sin(s), (T)0);
      return std::make_tuple(val, deriv);
    } else if (ver == ALGEBRAIC) {
      auto t = alpha() * (psi - psi_0()) * (real)M_1_PI, a = abs(t);
      auto mask = (a < (real)M_PI);
      auto x = (real)2 * a * (a - (real)1);
      auto val = select(mask, (real)1.0 - (a * a) / (x + (real)1), (T)0);
      auto dval_da = x / ((x + 1) * (x + 1));
      auto dval_dpsi = select(mask, alpha() * (real)M_1_PI * dval_da, (T)0);
      dval_dpsi = select(t < (real)0, -dval_dpsi, dval_dpsi);
      return std::make_tuple(val, dval_dpsi);
    } else {
      return std::make_tuple((T)0, (T)0);
    }
  }

  template <typename T>
  auto reg_call(T const &psi) const {
    if (version() == COSINE) {
      auto s = alpha() * (psi - psi_0());
      if (abs(s) >= M_PI) {
        return std::make_tuple((T)0, (T)0);
      } else {
        auto val = (real)0.5f * cos(s) + (real)0.5f;
        auto deriv = -(real)0.5f * alpha() * sin(s);
        return std::make_tuple(val, deriv);
      }
    } else if (version() == ALGEBRAIC) {
      auto t = alpha() * (psi - psi_0()) * (real)M_1_PI, a = abs(t);
      if (a >= M_PI) {
        return std::make_tuple((T)0, (T)0);
      } else {
        auto x = (real)2 * a * (a - (real)1);
        auto val = (real)1.0 - (a * a) / (x + (real)1);
        auto dval_da = x / ((x + (real)1) * (x + (real)1));
        auto dval_dpsi = alpha() * (real)M_1_PI * dval_da;
        if (t < (real)0)
          dval_dpsi = -dval_dpsi;
        return std::make_tuple(val, dval_dpsi);
      }
    } else {
      return std::make_tuple((T)0, (T)0);
    }
  }
};

class lambda : public lambda_expr<lambda> {
public:
  INST(lambda, FIELD(real, psi_0), FIELD(real, alpha),
       FIELD(lambda_version, version))

  lambda() : lambda((real)0, (real)0, COSINE){};
};
} // namespace cg::pid
