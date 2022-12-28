#pragma once
#include "../def/cuda_allocator.h"
#include "tuple.h"

namespace nitro::ind {
template <bool Indexed, typename T>
struct _cuda_allocator_impl;

template <typename T>
using cuda_allocator = typename _cuda_allocator_impl<is_indexed_v<T>, T>::type;

template <typename T>
struct _cuda_allocator_impl<false, T> {
  using type = def::cuda_allocator<T>;
};

// template <>
// struct _allocator_impl<false, bool> {
//   using type = bit::allocator;
// };

template <typename T>
struct _cuda_allocator_impl<true, T> {
  template <typename Idxes>
  struct _1;

  template <std::size_t... Idxes>
  struct _1<ind_seq<Idxes...>> {
    template <typename Types>
    struct _2;

    template <typename... Types>
    struct _2<type_list<Types...>> {
      class impl : public tuple<cuda_allocator<Types>...> {
      public:
        impl() : tuple<cuda_allocator<Types>...>{cuda_allocator<Types>()...} {}
      };
    };
  };

  using type = typename _1<idxes_t<T>>::template _2<subtypes_t<T>>::impl;
};
} // namespace nitro::ind