#pragma once
#include "def.h"
#include "par.h"

namespace nitro {
template <typename T, bool Indexed> struct const_iterator_aux;

template <typename T>
struct const_iterator_impl : const_iterator_aux<T, is_indexed<T>> {};

template <typename T> struct const_iterator_aux<T, false> {
  using type = def_const_iterator<T>;
};

template <typename T> struct const_iterator_aux<T, true> {
  template <size_t... ISeq>
  static auto aux(std::index_sequence<ISeq...>)
      -> par_const_iterator<T, ISeq...> {}

  using type = decltype(aux(std::make_index_sequence<T::num_types>{}));
};
} // namespace nitro