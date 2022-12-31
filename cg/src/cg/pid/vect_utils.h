#pragma once
#include <cg/pid/eval_forces.h>
#include <vcl/vectormath_trig.h>

namespace cg::pid {
// template <typename Mask, typename E, typename = void>
// struct select_impl {
//   static auto impl(Mask const &mask, E const &if_true, E const &if_false) {
//     return ::select(mask, if_true, if_false);
//   }
// };
//
// template <typename Mask, typename E>
// auto select_(Mask const &mask, E const &if_true, E const &if_false) {
//   return select_impl<Mask, E>::impl(mask, if_true, if_false);
// }
//
// template <typename Mask, typename E>
// struct select_impl<Mask, E, std::enable_if_t<vect::is_indexed_v<E>>> {
//   template <std::size_t... Idxes>
//   static auto aux(Mask const &mask, E const &if_true, E const &if_false,
//                   vect::ind_seq<Idxes...>) {
//     return E(select_(mask, if_true.template get<Idxes>(),
//                      if_false.template get<Idxes>())...);
//   }
//
//   static auto impl(Mask const &mask, E const &if_true, E const &if_false) {
//     return aux(mask, if_true, if_false, vect::idxes_t<E>{});
//   }
// };

template <typename T, std::size_t N, std::size_t W>
auto to_array(vect::lane<T, N, W> const &lane) {
  vect::array<T, N> elems;
  elems.template at_lane<N, W>(0) = lane;
  return elems;
}
} // namespace cg::pid