#pragma once
#include "../expr.h"
#include "../tuple.h"
#include "decl.h"
#include <cstddef>
#include <utility>

namespace nitro {

template <typename T, size_t N, size_t... ISeq>
class par_lane_const_at_expr
    : public auto_expr<par_lane_const_at_expr<T, N, ISeq...>, T>,
      public tuple_wrapper<
          lane_const_at_expr<typename T::template ith_type<ISeq>, N>...> {
public:
  using Base = tuple_wrapper<
      lane_const_at_expr<typename T::template ith_type<ISeq>, N>...>;
  using Base::Base;
  using Base::get;

  par_lane_const_at_expr(par_lane_const_at_expr const &other)
      : Base(static_cast<Base const &>(other)) {}

  par_lane_const_at_expr(par_lane_const_at_expr &&other) noexcept
      : Base(static_cast<Base &&>(other)) {}

  par_lane_const_at_expr &
  operator=(par_lane_const_at_expr const &other) = delete;
  par_lane_const_at_expr &operator=(par_lane_const_at_expr &&other) = delete;
};
} // namespace nitro