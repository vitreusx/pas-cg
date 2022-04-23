#pragma once
#include "def.h"
#include "par.h"

namespace nitro {

template <typename T, bool Indexed> struct allocator_aux;

template <typename T>
struct allocator_impl : allocator_aux<T, is_indexed<T>> {};

template <typename T> struct allocator_aux<T, false> {
  using type = def_allocator<T>;
};

template <typename T> struct allocator_aux<T, true> {
  template <size_t... ISeq>
  static auto aux(std::index_sequence<ISeq...>) -> par_allocator<T, ISeq...> {}

  using type = decltype(aux(std::make_index_sequence<T::num_types>{}));
};
} // namespace nitro