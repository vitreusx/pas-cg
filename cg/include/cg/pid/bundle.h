#pragma once
#include <cg/types/amp.h>

namespace cg::pid {
template <typename E> struct bundle_expr : public nitro::ind_expr<E> {
  EXPR_BODY(i1, i2, orig_dist, type)
};

template <typename E> struct bundle_auto_expr : public bundle_expr<E> {
  AUTO_EXPR_BODY(i1, i2, orig_dist, type)
};

using bundle_base = nitro::tuple_wrapper<int, int, real, int16_t>;

class bundle : public bundle_auto_expr<bundle>, public bundle_base {
public:
  using Base = bundle_base;
  using Base::Base;
  using Base::get;
};
} // namespace cg::pid

namespace nitro {
template <> struct is_indexed_impl<cg::pid::bundle> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::pid::bundle> {
  using type = cg::pid::bundle_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::pid::bundle> {
  using type = cg::pid::bundle_auto_expr<E>;
};
}