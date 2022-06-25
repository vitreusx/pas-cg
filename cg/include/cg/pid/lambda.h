#pragma once
#include <cg/types/amp.h>
#include <cg/vect/vect.h>

namespace cg::pid {
enum lambda_version {
  COSINE,
  ALGEBRAIC
};

template <typename E> struct lambda_expr {
  EXPR(psi_0, alpha, version)

  bool supp(real psi) const {
    return abs(alpha() * (psi - psi_0())) < M_PI;
  }

  inline std::tuple<real, real> operator()(real psi) const {
    switch (version()) {
    case COSINE: {
      auto s = alpha() * (psi - psi_0());
      if (abs(s) >= M_PI) {
        return std::make_tuple((real)0, (real)0);
      } else {
        auto val = 0.5f * cos(s) + 0.5f;
        auto deriv = -0.5f * alpha() * sin(s);
        return std::make_tuple(val, deriv);
      }
    }
    case ALGEBRAIC: {
      auto t = alpha() * (psi - psi_0()) * M_1_PI, a = abs(t);
      if (a >= M_PI) {
        return std::make_tuple((real)0, (real)0);
      } else {
        auto x = 2 * a * (a - 1);
        auto val = 1.0 - (a * a) / (x + 1);
        auto dval_da = x / ((x + 1) * (x + 1));
        auto dval_dpsi = alpha() * M_1_PI * dval_da;
        if (t < 0)
          dval_dpsi = -dval_dpsi;
        return std::make_tuple(val, dval_dpsi);
      }
    }
    default:
      return std::make_tuple((real)0, (real)0);
    }
  }
};

class lambda : public lambda_expr<lambda> {
public:
  INST(lambda, FIELD(real, psi_0), FIELD(real, alpha),
       FIELD(lambda_version, version))

  lambda() : lambda(0, 0, COSINE){};
};
} // namespace cg::pid
