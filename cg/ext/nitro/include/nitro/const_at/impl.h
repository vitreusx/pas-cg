#pragma once
#include "def.h"
#include "par.h"

namespace nitro {

template <typename T, bool Indexed> struct const_at_expr_aux;

template <typename T>
struct const_at_expr_impl : const_at_expr_aux<T, is_indexed<T>> {};

template <typename T> struct const_at_expr_aux<T, false> {
  using type = def_const_at<T>;
};

template <typename T> struct const_at_expr_aux<T, true> {
  template <size_t... ISeq>
  static auto aux(std::index_sequence<ISeq...>)
      -> par_const_at_expr<T, ISeq...> {}

  using type = decltype(aux(std::make_index_sequence<T::num_types>{}));
};
} // namespace nitro