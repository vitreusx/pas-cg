#pragma once
#include <cstddef>

namespace nitro {
template <typename... Types> struct type_list;

template <size_t I, typename Head, typename... Tail>
struct ith_type_impl : ith_type_impl<I - 1, Tail...> {};

template <typename Head, typename... Tail>
struct ith_type_impl<0, Head, Tail...> {
  using type = Head;
};

template <typename... Types> struct type_list {
  static constexpr size_t num_types = sizeof...(Types);

  template <size_t I>
  using ith_type = typename ith_type_impl<I, Types...>::type;
};

template <typename T, typename List> struct append_impl;

template <typename T, typename List>
using append = typename append_impl<T, List>::type;

template <typename T, typename... Rest>
struct append_impl<T, type_list<Rest...>> {
  using type = type_list<T, Rest...>;
};

template <typename List, typename T> struct prepend_impl;

template <typename List, typename T>
using prepend = typename prepend_impl<List, T>::type;

template <typename T, typename... Rest>
struct prepend_impl<type_list<Rest...>, T> {
  using type = type_list<Rest..., T>;
};

template <typename List1, typename List2> struct join_impl;

template <typename List1, typename List2>
using join = typename join_impl<List1, List2>::type;

template <typename List1> struct join_impl<List1, type_list<>> {
  using type = List1;
};

template <typename List1, typename Head2, typename... Tail2>
struct join_impl<List1, type_list<Head2, Tail2...>> {
  using type = join<prepend<List1, Head2>, type_list<Tail2...>>;
};

template <typename Types> struct concat_impl;

template <typename Types> using concat = typename concat_impl<Types>::type;

template <> struct concat_impl<type_list<>> { using type = type_list<>; };

template <typename Head, typename... Tail>
struct concat_impl<type_list<Head, Tail...>> {
  using type = join<Head, concat<type_list<Tail...>>>;
};

template <typename F, typename List> struct map_impl;

template <typename F, typename List>
using map = typename map_impl<F, List>::type;

template <typename F, typename... Types>
struct map_impl<F, type_list<Types...>> {
  using type = type_list<typename F::template apply<Types>...>;
};

template <template <typename> typename F> struct to_functor {
  template <typename T> using apply = F<T>;
};
}
