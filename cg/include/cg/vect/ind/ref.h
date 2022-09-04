#pragma once
#include "../bit/ref.h"
#include "ind_seq.h"
#include "tuple.h"
#include "type_list.h"
#include "type_traits.h"

namespace nitro::ind {
template <bool Indexed, typename T>
struct _ref_impl;

template <typename T>
using ref = typename _ref_impl<is_indexed_v<T>, T>::type;

template <typename T>
struct _ref_impl<false, T> {
  using type = T &;
};

// template <>
// struct _ref_impl<false, bool> {
//   using type = bit::ref;
// };

template <typename T>
struct _ref_impl<true, T> {
  template <typename Idxes>
  struct _1;

  template <std::size_t... Idxes_>
  struct _1<ind_seq<Idxes_...>> {
    template <typename Types>
    struct _2;

    template <typename... Types>
    struct _2<type_list<Types...>> {
      class impl : public expr_t<T, impl>, public tuple<ref<Types>...> {
      public:
        using tuple<ref<Types>...>::tuple;
        using tuple<ref<Types>...>::get;
        using tuple<ref<Types>...>::operator=;
        using expr_t<T, impl>::Idxes;

        impl(impl const &other) = default;

        auto &operator=(impl const &other) {
          tuple<ref<Types>...>::operator=(other);
          return *this;
        }

        impl(T &inst) : tuple<ref<Types>...>(inst.template get<Idxes_>()...) {}

        friend void swap(impl const &x, impl const &y) {
          using std::swap;
          (..., swap(x.template get<Idxes_>(), y.template get<Idxes_>()));
        }
      };
    };
  };

  using type = typename _1<idxes_t<T>>::template _2<subtypes_t<T>>::impl;
};
} // namespace nitro::ind