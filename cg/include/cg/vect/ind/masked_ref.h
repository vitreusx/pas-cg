#pragma once
// #include "../bit/lane_ref.h"
#include "../def/masked_ref.h"
#include "ind_seq.h"
#include "tuple.h"
#include "type_list.h"
#include "type_traits.h"

namespace nitro::ind {
template <bool Indexed, typename T, typename Idxes, typename Mask>
struct _masked_ref_impl;

template <typename T, typename Idxes, typename Mask>
using masked_ref =
    typename _masked_ref_impl<is_indexed_v<T>, T, Idxes, Mask>::type;

template <typename T, typename Idxes, typename Mask>
struct _masked_ref_impl<false, T, Idxes, Mask> {
  using type = def::masked_ref<T, Idxes, Mask>;
};

// template <std::size_t N, std::size_t W>
// struct _lane_ref_impl<false, bool, N, W> {
//   using type = bit::lane_ref<N, W>;
// };

template <typename T, typename Idx, typename Mask>
struct _masked_ref_impl<true, T, Idx, Mask> {
  template <typename Idxes>
  struct _1;

  template <std::size_t... _Idxes>
  struct _1<ind_seq<_Idxes...>> {
    template <typename Types>
    struct _2;

    template <typename... Types>
    struct _2<type_list<Types...>> {
      class impl : public expr_t<T, impl>,
                   public tuple<masked_ref<Types, Idx, Mask>...> {
      public:
        using Base = tuple<masked_ref<Types, Idx, Mask>...>;
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