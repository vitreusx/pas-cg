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
      using _Types = type_list<Types...>;

      template <std::size_t I>
      using _subtype = lane<std::tuple_element_t<I, _Types>, N, W>;

      class impl : public expr_t<T, impl>, public tuple<lane<Types, N, W>...> {
      public:
        using Base = tuple<lane<Types, N, W>...>;
        using Base::Base;
        using Base::get;
        using Base::operator=;
        using expr_t<T, impl>::Idxes;

        template <typename E, typename = std::enable_if_t<is_indexed_v<E>>>
        impl(E const &e)
            : Base{_subtype<_Idxes>(e.template get<_Idxes>())...} {}
      };
    };
  };

  using type = typename _1<idxes_t<T>>::template _2<subtypes_t<T>>::impl;
};

#ifdef __AVX512F___
inline constexpr int VECT_BITS = 512;
#else
inline constexpr int VECT_BITS = 256;
#endif
inline constexpr int VECT_BYTES = VECT_BITS / 8;

} // namespace nitro::ind