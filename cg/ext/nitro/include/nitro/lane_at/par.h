#pragma once
#include "../expr.h"
#include "../tuple.h"
#include "decl.h"
#include <cstddef>
#include <utility>

namespace nitro {

template <typename T, size_t N, size_t... ISeq>
class par_lane_at_expr
    : public auto_expr<par_lane_at_expr<T, N, ISeq...>, T>,
      public tuple_wrapper<
          lane_at_expr<typename T::template ith_type<ISeq>, N>...> {
public:
  using Base =
      tuple_wrapper<lane_at_expr<typename T::template ith_type<ISeq>, N>...>;
  using Base::Base;
  using Base::get;

  par_lane_at_expr(par_lane_at_expr const &other)
      : Base(static_cast<Base const &>(other)) {}

  par_lane_at_expr(par_lane_at_expr &&other) noexcept
      : Base(static_cast<Base &&>(other)) {}

  par_lane_at_expr &operator=(par_lane_at_expr const &other) {
    (..., assign<ISeq>(other));
    return *this;
  }

  template <typename E>
  par_lane_at_expr &operator=(ind_expr<E> const &e) {
    (..., assign<ISeq>(static_cast<E const &>(e)));
    return *this;
  }

  par_lane_at_expr &operator=(par_lane_at_expr &&other) noexcept {
    (..., assign<ISeq>(other));
    return *this;
  }

  template <typename E>
  par_lane_at_expr &operator=(ind_expr<E> &&e) {
    (..., assign<ISeq>(static_cast<E &&>(e)));
    return *this;
  }

private:
  template <size_t I, typename E>  void assign(E const &e) {
    this->template get<I>() = e.template get<I>();
  }
};
} // namespace nitro