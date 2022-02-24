#pragma once
#include <cg/types/amp.h>

namespace cg::vafm {
template <typename E> struct tip_expr : public nitro::ind_expr<E> {
  EXPR_BODY(res_idx, afm_orig, afm_vel)
};

template <typename E> struct tip_auto_expr : public tip_expr<E> {
  AUTO_EXPR_BODY(res_idx, afm_orig, afm_vel)
};

using tip_base = nitro::tuple_wrapper<int, vec3r, vec3r>;

class tip : public tip_expr<tip>, public tip_base {
public:
  using Base = tip_base;
  using Base::Base;
  using Base::get;
};
} // namespace cg::vafm

namespace nitro {
template <> struct is_indexed_impl<cg::vafm::tip> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::vafm::tip> {
  using type = cg::vafm::tip_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::vafm::tip> {
  using type = cg::vafm::tip_auto_expr<E>;
};
}; // namespace nitro