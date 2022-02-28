#pragma once
#include "contact_type.h"
#include <cg/amino/sync_data.h>
#include <cg/types/amp.h>

namespace cg::qa {
enum contact_status { FORMING_OR_FORMED, BREAKING };

template <typename E> struct contact_expr : public nitro::ind_expr<E> {
  EXPR_BODY(i1, i2, type, status, ref_time, sync_diff1, sync_diff2);
};

template <typename E> struct contact_auto_expr : public contact_expr<E> {
  AUTO_EXPR_BODY(i1, i2, type, status, ref_time, sync_diff1, sync_diff2);
};

using contact_base =
    nitro::tuple_wrapper<int, int, contact_type, contact_status, real,
                         sync_data, sync_data>;

class contact : public contact_auto_expr<contact>, public contact_base {
public:
  using Base = contact_base;
  using Base::Base;
  using Base::get;
};
} // namespace cg::qa

namespace nitro {
template <> struct is_indexed_impl<cg::qa::contact> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::qa::contact> {
  using type = cg::qa::contact_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::qa::contact> {
  using type = cg::qa::contact_auto_expr<E>;
};
}; // namespace nitro