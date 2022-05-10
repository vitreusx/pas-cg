#pragma once
#include "../aa_heur_pair.h"
#include <cg/types/amp.h>

namespace cg::heur_dih {
template <typename E> struct heur_dih_expr : public nitro::ind_expr<E> {
  EXPR_BODY(i1, i2, i3, i4, type)
};

template <typename E> struct heur_dih_auto_expr : public heur_dih_expr<E> {
  AUTO_EXPR_BODY(i1, i2, i3, i4, type)
};

using heur_dih_base = nitro::tuple_wrapper<int, int, int, int, aa_heur_pair>;

class heur_dih : public heur_dih_expr<heur_dih>, public heur_dih_base {
public:
  using Base = heur_dih_base;
  using Base::Base;
  using Base::get;
};
} // namespace cg::heur_dih

namespace nitro {
template <>
struct is_indexed_impl<cg::heur_dih::heur_dih> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::heur_dih::heur_dih> {
  using type = cg::heur_dih::heur_dih_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::heur_dih::heur_dih> {
  using type = cg::heur_dih::heur_dih_auto_expr<E>;
};
} // namespace nitro