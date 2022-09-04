#pragma once
#include <cstdint>
#include <tuple>

namespace nitro::ind {
template <typename... Types>
class type_list {};
} // namespace nitro::ind

namespace std {
template <size_t I, typename Head, typename... Tail>
struct tuple_element<I, nitro::ind::type_list<Head, Tail...>>
    : tuple_element<I - 1, nitro::ind::type_list<Tail...>> {};

template <typename Head, typename... Tail>
struct tuple_element<0, nitro::ind::type_list<Head, Tail...>> {
  using type = Head;
};
}