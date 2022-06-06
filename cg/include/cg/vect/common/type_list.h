#pragma once
#include <cstdint>
#include <tuple>

template <typename... Types>
class type_list {};

namespace std {
template <size_t I, typename Head, typename... Tail>
struct tuple_element<I, type_list<Head, Tail...>>
    : tuple_element<I - 1, type_list<Tail...>> {};

template <typename Head, typename... Tail>
struct tuple_element<0, type_list<Head, Tail...>> {
  using type = Head;
};
} // namespace std