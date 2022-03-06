#pragma once
#include "contact_type.h"
#include <cg/amino/sync_data.h>
#include <cg/types/amp.h>

namespace cg::qa {
template <typename E> struct candidate_expr : public nitro::ind_expr<E> {
  EXPR_BODY(i1, i2, dist, free_pair_idx, type, sync_diff1, sync_diff2)
};

template <typename E> struct candidate_auto_expr : public candidate_expr<E> {
  AUTO_EXPR_BODY(i1, i2, dist, free_pair_idx, type, sync_diff1, sync_diff2)
};

using candidate_base = nitro::tuple_wrapper<int, int, real, int, contact_type,
                                            sync_data, sync_data>;

class candidate : public candidate_auto_expr<candidate>, public candidate_base {
public:
  using Base = candidate_base;
  using Base::Base;
  using Base::get;
};
} // namespace cg::qa

namespace nitro {
template <>
struct is_indexed_impl<cg::qa::candidate> : public std::true_type {};

template <typename E> struct expr_impl<E, cg::qa::candidate> {
  using type = cg::qa::candidate_expr<E>;
};

template <typename E> struct auto_expr_impl<E, cg::qa::candidate> {
  using type = cg::qa::candidate_auto_expr<E>;
};
}; // namespace nitro