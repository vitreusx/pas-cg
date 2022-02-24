#pragma once
#include <cg/types/amp.h>
#include <nitro/nitro.h>

namespace cg::pid {
enum lambda_version { COSINE, ALGEBRAIC };

template <typename E> struct lambda_expr : public nitro::ind_expr<E> {
  EXPR_BODY(psi_0, alpha, version);

  bool supp(real psi) const { return abs(alpha() * (psi - psi_0())) < M_PI; }

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
    }
  }
};

template <typename E> struct lambda_auto_expr : public lambda_expr<E> {
  AUTO_EXPR_BODY(psi_0, alpha, version);
};

using lambda_base = nitro::tuple_wrapper<int, int, lambda_version>;

class lambda : public lambda_auto_expr<lambda>, public lambda_base {
public:
  using Base = lambda_base;
  using Base::Base;
  using Base::get;

  lambda() : Base(0, 0, COSINE){};
};
} // namespace cg::pid

namespace nitro {
template <> struct is_indexed_impl<cg::pid::lambda> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::pid::lambda> {
  using type = cg::pid::lambda_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::pid::lambda> {
  using type = cg::pid::lambda_auto_expr<E>;
};
}; // namespace nitro