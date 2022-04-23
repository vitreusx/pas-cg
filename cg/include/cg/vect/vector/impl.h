#pragma once
#include "decl.h"
#include "def.h"
#include "par.h"

namespace nitro {
template <typename T, typename Alloc, typename Idx, bool Indexed>
struct vector_impl_aux;

template <typename T, typename Alloc, typename Idx>
struct vector_impl : vector_impl_aux<T, Alloc, Idx, is_indexed<T>> {};

template <typename T, typename Alloc, typename Idx>
struct vector_impl_aux<T, Alloc, Idx, false> {
  using type = def_vector<T, Alloc, Idx>;
};

template <typename T, typename Alloc, typename Idx>
struct vector_impl_aux<T, Alloc, Idx, true> {
  template <size_t... ISeq>
  static auto aux(std::index_sequence<ISeq...>)
      -> par_vector<T, Alloc, Idx, ISeq...> {}

  using type = decltype(aux(std::make_index_sequence<T::num_types>{}));
};
} // namespace nitro
