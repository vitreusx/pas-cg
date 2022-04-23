#pragma once
#include "def.h"
#include "par.h"

namespace nitro {

template <typename T, size_t N, bool Indexed> struct lane_const_at_aux;

template <typename T, size_t N>
struct lane_const_at_expr_impl : lane_const_at_aux<T, N, is_indexed<T>> {};

template <typename T, size_t N> struct lane_const_at_aux<T, N, false> {
  using type = def_lane_const_at<T, N>;
};

template <typename T, size_t N> struct lane_const_at_aux<T, N, true> {
  template <size_t... ISeq>
  static auto aux(std::index_sequence<ISeq...>)
      -> par_lane_const_at_expr<T, ISeq...> {}

  using type = decltype(aux(std::make_index_sequence<T::num_types>{}));
};
} // namespace nitro