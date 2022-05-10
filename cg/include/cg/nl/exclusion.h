#pragma once
#include <cg/types/amp.h>

namespace cg::nl {
template <typename E> struct exclusion_expr : public nitro::ind_expr<E> {
  EXPR_BODY(i1, i2)
};

template <typename E> struct exclusion_auto_expr : public exclusion_expr<E> {
  AUTO_EXPR_BODY(i1, i2)
};

using exclusion_base = nitro::tuple_wrapper<int, int>;

class exclusion : public exclusion_auto_expr<exclusion>, public exclusion_base {
public:
  using Base = exclusion_base;
  using Base::Base;
  using Base::get;
};
} // namespace cg::nl

namespace nitro {
template <>
struct is_indexed_impl<cg::nl::exclusion> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::nl::exclusion> {
  using type = cg::nl::exclusion_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::nl::exclusion> {
  using type = cg::nl::exclusion_auto_expr<E>;
};
}