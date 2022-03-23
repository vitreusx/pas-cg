#pragma once
#include "../expr.h"
#include "../tuple.h"
#include "decl.h"
#include <cstddef>
#include <utility>

namespace nitro {

template <typename Types, size_t... ISeq>
class par_at_expr : public auto_expr<par_at_expr<Types, ISeq...>, Types>,
                    public tuple_wrapper<
                        at_expr<typename Types::template ith_type<ISeq>>...> {
public:
  using Base =
      tuple_wrapper<at_expr<typename Types::template ith_type<ISeq>>...>;
  using Base::Base;
  using Base::get;

   par_at_expr(par_at_expr const &other)
      : Base(static_cast<Base const &>(other)) { }

   par_at_expr(par_at_expr &&other)
      : Base(static_cast<Base &&>(other)) {}

   par_at_expr &operator=(par_at_expr const &other) {
    (..., assign<ISeq>(other));
    return *this;
  }

  template <typename E>
   par_at_expr &operator=(ind_expr<E> const &e) {
    (..., assign<ISeq>(static_cast<E const &>(e)));
    return *this;
  }

   par_at_expr &operator=(par_at_expr &&other) noexcept {
    (..., move<ISeq>(other));
    return *this;
  }

  template <typename E>
   par_at_expr &operator=(ind_expr<E> &&e) noexcept {
    (..., move<ISeq>(static_cast<E &&>(e)));
    return *this;
  }

  friend void swap(par_at_expr<Types, ISeq...> const& u, par_at_expr<Types, ISeq...> const& v) {
    using std::swap;
    (..., swap(u.template get<ISeq>(), v.template get<ISeq>()));
  }

private:
  template <size_t I, typename E>  void assign(E const &e) {
    this->template get<I>() = e.template get<I>();
  }

  template <size_t I, typename E>  void move(E &&e) noexcept {
    this->template get<I>() = std::move(e.template get<I>());
  }
};
} // namespace nitro