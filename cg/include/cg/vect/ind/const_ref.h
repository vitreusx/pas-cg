#pragma once
// #include "../bit/const_ref.h"
#include "ind_seq.h"
#include "tuple.h"
#include "type_list.h"
#include "type_traits.h"

namespace nitro::ind {
template <bool Indexed, typename T>
struct _const_ref_impl;

template <typename T>
using const_ref = typename _const_ref_impl<is_indexed_v<T>, T>::type;

template <typename T>
struct _const_ref_impl<false, T> {
  using type = T const &;
};

// template <>
// struct _const_ref_impl<false, bool> {
//   using type = bit::const_ref;
// };

template <typename T>
struct _const_ref_impl<true, T> {
  template <typename Idxes>
  struct _1;

  template <std::size_t... Idxes_>
  struct _1<ind_seq<Idxes_...>> {
    template <typename Types>
    struct _2;

    template <typename... Types>
    struct _2<type_list<Types...>> {
      class impl : public expr_t<T, impl>, public tuple<const_ref<Types>...> {
      public:
        using tuple<const_ref<Types>...>::tuple;
        using tuple<const_ref<Types>...>::get;
        using expr_t<T, impl>::Idxes;

        impl(T const &inst)
            : tuple<const_ref<Types>...>(inst.template get<Idxes_>()...) {}
      };
    };
  };

  using type = typename _1<idxes_t<T>>::template _2<subtypes_t<T>>::impl;
};
} // namespace nitro::ind