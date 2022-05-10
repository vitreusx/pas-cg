#pragma once
#include "../aa_heur_pair.h"
#include <cg/types/amp.h>

namespace cg::heur_ang {
template <typename E> struct heur_ang_expr : public nitro::ind_expr<E> {
  EXPR_BODY(i1, i2, i3, type)
};

template <typename E> struct heur_ang_auto_expr : public heur_ang_expr<E> {
  AUTO_EXPR_BODY(i1, i2, i3, type)
};

using heur_ang_base = nitro::tuple_wrapper<int, int, int, aa_heur_pair>;

class heur_ang : public heur_ang_expr<heur_ang>, public heur_ang_base {
public:
  using Base = heur_ang_base;
  using Base::Base;
  using Base::get;
};
} // namespace cg::heur_ang

namespace nitro {
template <>
struct is_indexed_impl<cg::heur_ang::heur_ang> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::heur_ang::heur_ang> {
  using type = cg::heur_ang::heur_ang_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::heur_ang::heur_ang> {
  using type = cg::heur_ang::heur_ang_auto_expr<E>;
};
} // namespace nitro