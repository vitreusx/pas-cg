#pragma once
// #include "../bit/lane_const_ref.h"
#include "../def/lane.h"
#include "ind_seq.h"
#include "tuple.h"
#include "type_list.h"
#include "type_traits.h"

namespace nitro::ind {
template <bool Indexed, typename T, std::size_t N, std::size_t W>
struct _lane_impl;

template <typename T, std::size_t N, std::size_t W>
using lane = typename _lane_impl<is_indexed_v<T>, T, N, W>::type;

template <typename T, std::size_t N, std::size_t W>
struct _lane_impl<false, T, N, W> {
  using type = def::lane<T, N, W>;
};

// template <std::size_t N, std::size_t W>
// struct _lane_const_ref_impl<false, bool, N, W> {
//   using type = bit::lane_const_ref<N, W>;
// };

template <typename T, std::size_t N, std::size_t W>
struct _lane_impl<true, T, N, W> {
  template <typename Idxes>
  struct _1;

  template <std::size_t... _Idxes>
  struct _1<ind_seq<_Idxes...>> {
    template <typename Types>
    struct _2;

    template <typename... Types>
    struct _2<type_list<Types...>> {
      class impl : public expr_t<T, impl>, public tuple<lane<Types, N, W>...> {
      public:
        using tuple<lane<Types, N, W>...>::tuple;
        using tuple<lane<Types, N, W>...>::get;
        using tuple<lane<Types, N, W>...>::operator=;
        using expr_t<T, impl>::Idxes;
      };
    };
  };

  using type = typename _1<idxes_t<T>>::template _2<subtypes_t<T>>::impl;
};

} // namespace nitro::ind