#pragma once
// #include "../bit/lane_ref.h"
// #include "../def/sparse_ref.h"
#include "../def/const_ref.h"
#include "ind_seq.h"
#include "tuple.h"
#include "type_list.h"
#include "type_traits.h"

namespace nitro::ind {
template <bool Indexed, typename T, typename Idxes>
struct _sparse_const_ref_impl;

template <typename T, typename Idxes>
using sparse_const_ref =
    typename _sparse_const_ref_impl<is_indexed_v<T>, T, Idxes>::type;

template <typename T, typename Idxes>
struct _sparse_const_ref_impl<false, T, Idxes> {
  using type = def::sparse_const_ref<T, Idxes>;
};

// template <std::size_t N, std::size_t W>
// struct _lane_ref_impl<false, bool, N, W> {
//   using type = bit::lane_ref<N, W>;
// };

template <typename T, typename Idx>
struct _sparse_const_ref_impl<true, T, Idx> {
  template <typename Idxes>
  struct _1;

  template <std::size_t... _Idxes>
  struct _1<ind_seq<_Idxes...>> {
    template <typename Types>
    struct _2;

    template <typename... Types>
    struct _2<type_list<Types...>> {
      class impl : public expr_t<T, impl>,
                   public tuple<sparse_const_ref<Types, Idx>...> {
      public:
        using Base = tuple<sparse_const_ref<Types, Idx>...>;
        using Base::Base;
        using Base::get;
        using Base::operator=;
        using expr_t<T, impl>::Idxes;
      };
    };
  };

  using type = typename _1<idxes_t<T>>::template _2<subtypes_t<T>>::impl;
};
} // namespace nitro::ind