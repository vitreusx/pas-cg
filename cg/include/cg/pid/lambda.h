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
      auto val = 0.5f * cos(s) + 0.5f;
      auto deriv = -0.5f * alpha() * sin(s);
      return std::make_tuple(val, deriv);
    }
    case ALGEBRAIC: {
      auto s = alpha() * (psi - psi_0());
      auto t = abs(s / M_PI);
      auto x_inv = 1.0f / (2.0f * t * t - 2.0f * t - 1);
      auto val = (t * t - 2.0f * t + 1.0f) * x_inv;
      auto deriv = (2.0f * t * (t - 1.0f)) * x_inv * x_inv / M_PI;
      deriv *= (s < 0.0f ? -1.0f : 1.0f);
      return std::make_tuple(val, deriv);
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
